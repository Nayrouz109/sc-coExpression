#corOnly_parallelV3_disease
# this script was modified from corOnly_parallelV0_disease

# the difference is that it implements Queue multiprocessing shared memory 
# this is to minimize memory usage 

import numpy as np
import scipy.sparse
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from joblib import Parallel, delayed
from pathlib import Path
import multiprocessing as mp




def rank(data):
    orig_shape = data.shape
    data = np.argsort(np.argsort(data, axis=None))
    return (data / (data.size - 1)).reshape(orig_shape)



def sparse_corr(A):
    """Compute Pearson correlation for a sparse matrix.
    
    Arguments:
        A {scipy.sparse matrix} -- Sparse matrix of gene expression values
    Returns:
        np.matrix -- Correlation matrix
    """
    N = A.shape[0]
    #C = ((A.T * A - (A.sum(axis=0).T * A.sum(axis=0) / N)) / (N - 1)).todense()
    C = ((A.T*A -(sum(A).T*sum(A)/N))/(N-1)).todense()
    V = np.sqrt(np.mat(np.diag(C)).T * np.mat(np.diag(C)))
    COR = np.divide(C, V + 1e-119)
    return COR

def process_patient(patient_id, 
                    adata, 
                    min_cells_threshold, 
                    correlation_type, 
                    replace_nans, 
                    queue):
    
    subset_data = adata[adata.obs['patientID'] == patient_id].X
    subset_data = scipy.sparse.csr_matrix(subset_data)

    if subset_data.shape[0] < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.")
        return

    if correlation_type == 'pearson':
        corr_matrix = sparse_corr(subset_data)
    elif correlation_type == 'spearman':
        corr_matrix, _ = spearmanr(subset_data, axis=0, nan_policy='omit')
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")

    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix)

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = np.nanmean(rank_normalized_matrix)

    # puts the rank-normalized matrix into the queue 
    queue.put(rank_normalized_matrix)
    print(f"Finished processing patient {patient_id}")



def aggregate_results(queue, 
                      num_results, 
                      gene_count):
    agg_matrix = csr_matrix((gene_count, gene_count), dtype=np.float64)
    count = 0 #keep track of how many matrices have been aggregated
    #num_results: total number of matrices expected to be retrieved & aggregated from the queue 

    while count < num_results:
        rank_matrix = queue.get()
        if rank_matrix is not None:
            agg_matrix += csr_matrix(rank_matrix)
            count += 1

    if count > 0:
        agg_matrix /= count

    return agg_matrix



def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             max_workers=6, 
                             analysis_type=None, 
                             categories=None, 
                             file_path=""):
    print("Starting correlation computation...")
    file_path = Path(file_path)
    queue = mp.Queue()
    print("hi")

    if analysis_type == "disease" and categories:
        subsets = {category: adata[adata.obs['disease'] == category] for category in categories}
        results = {}

        for category, subset in subsets.items():
            print(f"Processing subset for category: {category}")
            patient_ids = subset.obs['patientID'].unique()
            total_patients = len(patient_ids)
            print(f"Total number of patients in {category}: {total_patients}")

            processes = []
            for patient_id in patient_ids:
                p = mp.Process(target=process_patient, 
                               args=(patient_id, subset, min_cells_threshold, 
                                     correlation_type, replace_nans, queue))
                processes.append(p)
                p.start()

            for p in processes:
                p.join()

            aggregated_matrix = aggregate_results(queue, total_patients, subset.shape[1])
            summary_df = pd.DataFrame({
                "category": [category],
                "total_genes": [subset.shape[1]],
                "total_single_cells": [subset.shape[0]], 
                "total_patients": total_patients
            })

            results[category] = {
                'aggregated_rank_normalized_matrix': aggregated_matrix,
                'gene_names': subset.var_names.tolist(),
                'data_summary': summary_df,
                'file_path': str(file_path).replace(".h5ad", f"_{category}.h5ad")
            }

        return results

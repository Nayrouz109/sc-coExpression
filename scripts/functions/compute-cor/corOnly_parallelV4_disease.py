

### this is the same as V5 
### exception is how average agg was calculated 

import numpy as np
import scipy.sparse
import pandas as pd 
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from joblib import Parallel, delayed
import bottleneck 
from pathlib import Path
import bottleneck as bn

def rank(data):
    """Rank normalize data with NaN handling"""
    orig_shape = data.shape
    data = bn.nanrankdata(data).reshape(-1) - 1  # Ranks, ignoring NaNs
    # Normalize between 0 and 1
    return (data / np.nansum(~np.isnan(data))).reshape(orig_shape)

def sparse_corr(A):
    """Compute Pearson correlation for a sparse matrix.
    
    Arguments:
        A {scipy.sparse matrix} -- Sparse matrix of gene expression values
    Returns:
        np.matrix -- Correlation matrix
    """
    dense_matrix = A.toarray()
    COR = np.corrcoef(dense_matrix, rowvar=False)
    return COR

def process_patient(patient_id, adata, min_cells_threshold, correlation_type, replace_nans):
    """
    Process data for a single patient: compute correlation and rank normalize.

    Returns:
        tuple -- Patient ID and a dictionary with processed edge lists and rank-normalized matrix.
    """
    subset_data = adata[adata.obs['patientID'] == patient_id].X
    subset_data = scipy.sparse.csr_matrix(subset_data)

    if subset_data.shape[0] < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.")
        return patient_id, None
    
    # Compute correlation matrix
    if correlation_type == 'pearson':
        corr_matrix = sparse_corr(subset_data)
    elif correlation_type == 'spearman':
        corr_matrix, _ = spearmanr(subset_data, axis=0, nan_policy='omit')
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")
    
    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix)

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bottleneck.nanmean(rank_normalized_matrix)

    print(f"Finished processing patient {patient_id}")
    return patient_id, {
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }

def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             analysis_type=None, 
                             categories=None, 
                             file_path=""):
    print("Starting correlation computation...")
    file_path = Path(file_path)

    if analysis_type == "disease" and categories:
        subsets = {category: adata[adata.obs['disease'] == category] for category in categories}

        results = {}
        for category, subset in subsets.items():
            print(f"Processing subset for category: {category}")

            patient_ids = subset.obs['patientID'].unique()
            total_patients = len(patient_ids)
            print(f"Total number of patients in {category}: {total_patients}")

            aggregated_matrix = None
            gene_names = None

            for patient_id in patient_ids:
                print(f"Processing patient ID: {patient_id}")

                result = process_patient(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)

                if result is not None:
                    #patient_matrix = result['rank_normalized_matrix']
                    patient_matrix = result[1]["rank_normalized_matrix"]
                    if gene_names is None:
                        #gene_names = result['gene_names']
                        gene_names = result[1]["gene_names"]

                    # Aggregate the patient matrix into the aggregated matrix
                    if aggregated_matrix is None:
                        aggregated_matrix = patient_matrix
                    else:
                        aggregated_matrix = np.mean([aggregated_matrix, patient_matrix], axis=0)

                    print(f"Patient {patient_id} matrix aggregated.")
                else:
                    print(f"No valid matrix for patient {patient_id}.")

            # Convert aggregated matrix to a DataFrame with gene names as row and column names
            if aggregated_matrix is not None and gene_names is not None:
                aggregated_matrix = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

            summary_df = pd.DataFrame({
                "category": [category],
                "total_genes": [len(gene_names) if gene_names else 0],
                "total_single_cells": [subset.shape[0]], 
                "total_patients": total_patients
            })

            results[category] = {
                'aggregated_rank_normalized_matrix': aggregated_matrix,
                'gene_names': gene_names,
                'data_summary': summary_df,
                'file_path': str(file_path).replace(".h5ad", f"_{category}.h5ad")
            }

        return results

    else:
        # Non-disease analysis logic
        pass

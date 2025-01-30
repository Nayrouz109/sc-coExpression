


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
import multiprocessing





def matrix_to_edge_list_efficient(matrix, gene_names):
    """
    Convert a square matrix to an edge list DataFrame.

    Args:
        matrix (np.ndarray or np.matrix): Correlation or rank-normalized matrix.
        gene_names (list): List of gene names corresponding to matrix rows/columns.

    Returns:
        pd.DataFrame: Edge list containing gene1, gene2, and correlation/rank values.
    """
    # Ensure matrix is a numpy array
    if isinstance(matrix, np.matrix):
        matrix = np.asarray(matrix)

    # Extract upper triangular indices (excluding diagonal)
    non_zero_indices = np.triu_indices_from(matrix, k=1)
    edge_values = matrix[non_zero_indices]
    gene1 = [gene_names[i] for i in non_zero_indices[0]]
    gene2 = [gene_names[j] for j in non_zero_indices[1]]

    # Create the edge list DataFrame
    edge_list = pd.DataFrame({
        "gene1": gene1,
        "gene2": gene2,
        "value": edge_values.ravel()
    })

    # Filter out rows with 0 or NaN values
    edge_list = edge_list[(edge_list["value"] != 0) & edge_list["value"].notna()]

    return edge_list

def extract_metadata_from_filename(filename):
    """Extracts the study name from the filename."""
    name_parts = filename.split('_')
    cell_type = name_parts[0] 
    study_name = name_parts[-1].split('.')[0] 
    return cell_type, study_name


def rank(data):
    """Rank normalize data with NaN handling"""
    orig_shape = data.shape
    data = bn.nanrankdata(data).reshape(-1) - 1  # Ranks, ignoring NaNs
    # Normalize between 0 and 1
    return (data / np.nansum(~np.isnan(data))).reshape(orig_shape)

def rank_noNA(data):
    """Rank normalize data
    
    Rank standardize data to make nonparametric
    
    Arguments:
        data {np.array} -- 2-D coexpression network
    
    Returns:
        np.array -- Rank normalized between 0 and 1 array
    """
    orig_shape = data.shape
    #data = bottleneck.nanrankdata(data) - 1
    #return (data / np.sum(~np.isnan(data))).reshape(orig_shape)
    data = np.argsort(np.argsort(data, axis=None))  # Rank data
    return (data / (data.size - 1)).reshape(orig_shape)


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

def sparse_corrM(A):
    """manual Pearson correlation for a sparse matrix.
    
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
        #'corrank_edge_list': edgeList,
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }




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

    if analysis_type == "disease" and categories:
        subsets = {category: adata[adata.obs['disease'] == category] for category in categories}

        results = {}
        for category, subset in subsets.items():
            print(f"Processing subset for category: {category}")

            patient_ids = subset.obs['patientID'].unique()
            total_patients = len(patient_ids)
            print(f"total number of patients in {category}: {total_patients}")

            subset_results = Parallel(n_jobs=max_workers)(
                delayed(process_patient)(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)
                for patient_id in patient_ids
            )

            # this is to create a dictionary of numpy.matrices 
            rank_normalized_per_patientID = {}
            gene_names = None
            for patient_id, result in subset_results:
                if result is not None:
                    rank_normalized_per_patientID[patient_id] = result['rank_normalized_matrix']
                    if gene_names is None:
                        gene_names = result['gene_names']

            if rank_normalized_per_patientID:
                aggregated_matrix = np.mean(list(rank_normalized_per_patientID.values()), axis=0)
            else:
                aggregated_matrix = None

            aggregated_matrix = csr_matrix(aggregated_matrix)
            summary_df = pd.DataFrame({
                "category": [category],
                "total_genes": [len(gene_names)],
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
        # non-disease analysis logic
        pass



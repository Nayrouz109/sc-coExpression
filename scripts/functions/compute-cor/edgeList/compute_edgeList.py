

import numpy as np
import scipy.sparse
import pandas as pd 
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from joblib import Parallel, delayed
import bottleneck 
from pathlib import Path
from joblib import parallel_backend
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

    correlation_edge_list = matrix_to_edge_list_efficient(corr_matrix, adata.var_names.tolist())
    corrank_edge_list = matrix_to_edge_list_efficient(rank_normalized_matrix, adata.var_names.tolist())
    del corr_matrix
    del rank_normalized_matrix 

    # Merge edge lists
    edgeList = pd.merge(
        correlation_edge_list,
        corrank_edge_list,
        on=["gene1", "gene2"],
        how="inner",
        suffixes=("_cor", "_rank")
    )

    # Cleanup
    del correlation_edge_list
    del corrank_edge_list
    del subset_data 

    print(f"Finished processing patient {patient_id}")
    return patient_id, {'corrank_edge_list': edgeList}



def compute__edgeListCor(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             max_workers=7, 
                             analysis_type=None, 
                             file_path=""):
    """
    Compute correlation coefficients, rank normalize them, and return processed results.
    """
    print("Starting correlation computation for all patients...")
    file_path = Path(file_path)

    #cell_type, study_name = extract_metadata_from_filename(file_path.name)
    if analysis_type == "baseline":
        adata = adata[adata.obs['diagnosis'] == 'no']
    #cells = adata.shape[0]
    patient_ids = adata.obs['patientID'].unique()

    with parallel_backend('loky', n_jobs=max_workers):
        results = Parallel(n_jobs=max_workers)(
            delayed(process_patient)(patient_id, adata, min_cells_threshold, correlation_type, replace_nans)
            for patient_id in patient_ids
    )
    #results = Parallel(n_jobs=max_workers)(
        #delayed(process_patient)(patient_id, adata, min_cells_threshold, correlation_type, replace_nans)
        #for patient_id in patient_ids
    #)

    del adata
    rnkEdgeList_per_patientID = {}

    for patient_id, result in results:
        if result is not None:
            rnkEdgeList_per_patientID[patient_id] = result['corrank_edge_list']

    del results
    print("merging edge lists across patients")
    merged_df = pd.concat(
        {key: df.set_index(['gene1', 'gene2']) for key, df in rnkEdgeList_per_patientID.items()}, 
                          names=['patient', 'gene1', 'gene2']
                          ).reset_index()
    
    merged_df.set_index(['patient', 'gene1', 'gene2'], inplace=True)
    del rnkEdgeList_per_patientID

    print("Finished correlation computation.")


    return {
        'corEdgeList_per_patientID': merged_df,
    }

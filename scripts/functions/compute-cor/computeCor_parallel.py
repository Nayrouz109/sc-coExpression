

import numpy as np
import scipy.sparse
import pandas as pd 
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from joblib import Parallel, delayed
import bottleneck 
from pathlib import Path
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
    
    Arguments:
        patient_id {str} -- Patient ID
        adata {AnnData} -- AnnData object
        min_cells_threshold {int} -- Minimum number of cells required
        correlation_type {str} -- Type of correlation ('pearson' or 'spearman')
        replace_nans {bool} -- Whether to replace NaN values in the correlation matrix
    
    Returns:
        tuple -- Patient ID and a dictionary with correlation matrix and rank-normalized matrix
    """

    subset_data = adata[adata.obs['patientID'] == patient_id].X
    # Ensure sparse format is subscriptable
    subset_data = scipy.sparse.csr_matrix(subset_data)

    if subset_data.shape[0] < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.")
        return patient_id, None  # Skip patient
    
    # Step 1: Compute the correlation matrix based on correlation_type
    if correlation_type == 'pearson':
        print(f"\nComputing Pearson correlation for patient {patient_id}...")
        corr_matrix = sparse_corr(subset_data)  # Use sparse_corr for Pearson correlation
    elif correlation_type == 'spearman':
        corr_matrix, _ = spearmanr(subset_data, axis=0, nan_policy='omit')  # Compute Spearman correlation
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")
    
    # Step 2: Fill diagonal with 1 and rank normalize
    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix)

    # Step 3: Replace NaNs if needed
    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bottleneck.nanmean(rank_normalized_matrix)

    # Step 4: Convert matrices to edge lists
    print(f"Converting correlation matrix to edge list for patient {patient_id}...")
    correlation_edge_list = matrix_to_edge_list_efficient(corr_matrix, adata.var_names.tolist())

    # Step 4: Convert matrices to edge lists
    print(f"Converting rank correlation matrix to edge list for patient {patient_id}...")
    corrank_edge_list = matrix_to_edge_list_efficient(rank_normalized_matrix, adata.var_names.tolist())

    # Step 4: Return the results
    print(f"Finished processing patient {patient_id}")
    return patient_id, {'corrank_edge_list': corrank_edge_list, 
                        'correlation_edge_list': correlation_edge_list,
                        'rank_normalized_matrix': rank_normalized_matrix,
                        'gene_names': adata.var_names.tolist()} 




#def compute_gene_correlation(adata, min_cells_threshold=20):
def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             max_workers=7, 
                             analysis_type=None, 
                             file_path = ""):
    """
    Compute the correlation coefficients, rank normalize them,
    and return the results, leveraging parallelization.
    
    Arguments:
        adata {AnnData} -- The input AnnData object containing gene expression data and metadata
        correlation_type {str} -- Type of correlation ('pearson', 'spearman')
        replace_nans {bool} -- Replace NaNs in the correlation matrix with the median value
        min_cells_threshold {int} -- Minimum number of cells to compute correlation for each patient
        max_workers {int} -- Maximum number of workers for parallel computation
        analysis_type {str or None} -- Type of analysis. If 'baseline', subset `adata` where diagnosis == 'no`
    
    
    Returns:
        dict -- Dictionary containing:
            - 'correlation_per_patientID' : Raw correlation values per patientID
            - 'rank_normalized_per_patientID' : Rank normalized correlation matrix per patientID
            - 'aggregated_rank_normalized_matrix' : Aggregated rank normalized correlation matrix
            - 'gene_names': Gene names used for computation
    """
    print("SSStarting correlation computation for all patients...")
    file_path = Path(file_path)  # Ensure file_path is a Path object
    try:
        cell_type, study_name = extract_metadata_from_filename(file_path.name)
    except Exception as e:
        raise ValueError(f"Error extracting metadata from file path '{file_path}': {e}")
    

    # Subset AnnData if analysis_type == "baseline"
    if analysis_type == "baseline":
        print("Subsetting AnnData for baseline analysis: diagnosis == 'no'")
        if "diagnosis" not in adata.obs:
            raise ValueError("The 'diagnosis' column is missing in `adata.obs`. Cannot perform baseline analysis.")
        
        adata = adata[adata.obs['diagnosis'] == 'no']
        if adata.shape[0] == 0:
            raise ValueError("No cells available after subsetting for baseline analysis (diagnosis == 'no').")
        
    patient_ids = adata.obs['patientID'].unique()
    
    print(f"Using {max_workers} workers for parallel computation.")
    results = Parallel(n_jobs=max_workers)(
        delayed(process_patient)(patient_id, adata, min_cells_threshold, correlation_type, replace_nans)
        for patient_id in patient_ids
    )

    # Collect results
    #correlation_per_patientID = {}
    rnkEdgeList_per_patientID = {}
    corEdgeList_per_patientID = {}
    rank_normalized_per_patientID = {}
    gene_names = None  # Gene names should be the same across all patients, we take from the first patient

    for patient_id, result in results:
        if result is not None:
            #correlation_per_patientID[patient_id] = result['correlation_matrix']
            corEdgeList_per_patientID[patient_id] = result['correlation_edge_list']
            rnkEdgeList_per_patientID[patient_id] = result['corrank_edge_list']
            rank_normalized_per_patientID[patient_id] = result['rank_normalized_matrix']
            if gene_names is None:
                gene_names = result['gene_names']  # Store the gene names from the first patient
    

    corEdheList = pd.concat({key: df.set_index(['gene1', 'gene2']) 
                         for key, df in corEdgeList_per_patientID.items()}, 
                        names=['patient', 'gene1', 'gene2'])
    
    corRnkEdheList = pd.concat({key: df.set_index(['gene1', 'gene2']) 
                         for key, df in rnkEdgeList_per_patientID.items()}, 
                        names=['patient', 'gene1', 'gene2'])
    
    corEdgeList_per_patientID = corEdheList.reset_index()
    rnkEdgeList_per_patientID = corRnkEdheList.reset_index()

    del corEdheList
    del corRnkEdheList

    # Merge the two dataframes on 'patient', 'gene1', 'gene2'
    print("merging edge lists")
    merged_df = pd.merge(
        corEdgeList_per_patientID,
        rnkEdgeList_per_patientID,
        on=["patient", "gene1", "gene2"],
        how='inner', 
        suffixes=('_cor', '_rank')
    )
    merged_df.set_index(["patient", "gene1", "gene2"], inplace=True)
    #filtered_df = merged_df[merged_df['value_cor'] != 0]

    del corEdgeList_per_patientID
    del rnkEdgeList_per_patientID

    print(f"Completed processing for {len(rank_normalized_per_patientID)} patients.")

    # Aggregate the rank-normalized matrices by averaging across patients
    if rank_normalized_per_patientID:
        print("Aggregating rank-normalized matrices across all patients...")
        aggregated_matrix = np.mean(list(rank_normalized_per_patientID.values()), axis=0)
    else:
        aggregated_matrix = None

    aggregated_matrix = csr_matrix(aggregated_matrix)
    del rank_normalized_per_patientID

    summary_df = pd.DataFrame({
        "study_name": [study_name],
        "cell_type": [cell_type],
        "total_genes": [adata.shape[1]],
        "total_single_cells": [adata.shape[0]],
        "total_patient_IDs": [corEdheList.index.get_level_values('patient').nunique()]
    })

    print("Finished correlation computation.")
    del adata 
    return {
        #'correlation_per_patientID': correlation_per_patientID,
        'corEdgeList' :merged_df, 
        #'rank_normalized_per_patientID': rank_normalized_per_patientID,
        #'rnkEdgeList_per_patientID' : rnkEdgeList_per_patientID , 
        'aggregated_rank_normalized_matrix': aggregated_matrix, 
        'gene_names': gene_names,   # Return gene names along with the matrices
        'data_summary': summary_df
    }





import numpy as np
import scipy.sparse
import pandas as pd 
from scipy.stats import spearmanr
import bottleneck 
from pathlib import Path
import bottleneck as bn
import sys
import gc 
import psutil
import os 

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from saveCorDict import save_to_hdf5

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
    Process data for a single patient: compute correlation, rank normalize, and return weight.
    
    Returns:
        tuple -- Patient ID and a dictionary with processed edge lists and rank-normalized matrix,
                 along with the patient's weight (total single cells).
    """
    subset_adata = adata[adata.obs['patientID'] == patient_id]
    patient_cells = subset_adata.shape[0]  # Number of single cells in this patient
    
    if patient_cells < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.")
        return patient_id, None, 0  # Weight is zero for excluded patients
    
    # Compute correlation matrix
    if correlation_type == 'pearson':
        corr_matrix = sparse_corr(scipy.sparse.csr_matrix(subset_adata.X))
    elif correlation_type == 'spearman':
        corr_matrix, _ = spearmanr(subset_adata.X, axis=0, nan_policy='omit')
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")

    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix)

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bottleneck.nanmean(rank_normalized_matrix)

    print(f"Finished processing patient {patient_id}")

    del corr_matrix  
    return patient_id, {
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }, patient_cells  # Return patient weight


def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             category=None, 
                             file_path="", 
                             output_dir=""):
    print("Starting correlation computation...")
    file_path = Path(file_path)
    output_dir = Path(output_dir)

    print(f"Processing subset for category: {category}")

    subset = adata[adata.obs['disease'] == category]
    patient_ids = subset.obs['patientID'].unique()
    total_patients = len(patient_ids)
    print(f"Total number of patients in {category}: {total_patients}")

    aggregated_matrix = None
    gene_names = None
    patient_count = 0
    total_weight = 0  # Track total sum of patient weights

    for patient_id in patient_ids:
        print(f"Processing patient ID: {patient_id}")

        patient_id, result, patient_weight = process_patient(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)

        if result is not None:
            patient_matrix = result["rank_normalized_matrix"]
            if gene_names is None:
                gene_names = result["gene_names"]

            # Apply weight and aggregate
            if aggregated_matrix is None:
                aggregated_matrix = patient_matrix * patient_weight  # Weighted contribution
            else:
                aggregated_matrix += patient_matrix * patient_weight  # Weighted sum

            patient_count += 1
            total_weight += patient_weight  # Accumulate total weight
            print(f"Patient {patient_id} matrix aggregated with weight {patient_weight}.")
            del patient_matrix
            gc.collect()
        else:
            print(f"No valid matrix for patient {patient_id}.")
    
    # Compute weighted average
    if aggregated_matrix is not None and total_weight > 0:
        aggregated_matrix = aggregated_matrix / total_weight

    # Convert aggregated matrix to DataFrame
    if aggregated_matrix is not None and gene_names is not None:
        aggregated_matrixEdge = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

    summary_df = pd.DataFrame({
        "category": [category],
        "total_genes": [len(gene_names) if gene_names else 0],
        "total_single_cells": [total_weight],  
        "total_patients": patient_count
    })

    category_data = {
        'aggregated_rank_normalized_matrix': aggregated_matrix,
        'edgeList': aggregated_matrixEdge, 
        'gene_names': gene_names,
        'data_summary': summary_df, 
        'file_path': str(file_path).replace(".h5ad", f"_{category}.h5ad")
    }
    category_file_path = Path(category_data["file_path"])

    # Log memory usage
    process = psutil.Process(os.getpid())
    print(f"Memory usage before saving: {process.memory_info().rss / 1024 ** 2} MB")

    print(f"Saving results for category: {category}")
    save_to_hdf5(category_data, output_dir, category_file_path)
    print(f"Results for {category} saved successfully to {category_file_path}.")

    # Clean up memory
    del subset
    del aggregated_matrix
    del aggregated_matrixEdge
    del category_data
    del summary_df
    gc.collect()  # Final garbage collection
    print(f"Memory usage after cleanup: {process.memory_info().rss / 1024 ** 2} MB")

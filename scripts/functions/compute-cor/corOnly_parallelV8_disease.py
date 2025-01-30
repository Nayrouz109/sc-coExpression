import numpy as np
import scipy.sparse
import pandas as pd
from scipy.stats import spearmanr
import bottleneck as bn
from pathlib import Path
import sys
import gc
import psutil
import os
import h5py
import tracemalloc  # Import tracemalloc



def rank(data):
    """Rank normalize data with NaN handling"""
    orig_shape = data.shape
    data = bn.nanrankdata(data).reshape(-1) - 1  # Ranks, ignoring NaNs
    return (data / np.nansum(~np.isnan(data))).reshape(orig_shape)

def sparse_corr(A):
    """Compute Pearson correlation for a sparse matrix."""
    dense_matrix = A.toarray()
    COR = np.corrcoef(dense_matrix, rowvar=False)
    return COR

def process_patient(patient_id, adata, min_cells_threshold, correlation_type, replace_nans):
    """
    Process data for a single patient: compute correlation and rank normalize.
    """
    if adata[adata.obs['patientID'] == patient_id].X.shape[0] < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.")
        return patient_id, None
    
    tracemalloc.start_snapshot(f"Before processing patient {patient_id}")

    # Compute correlation matrix
    if correlation_type == 'pearson':
        corr_matrix = sparse_corr(scipy.sparse.csr_matrix(adata[adata.obs['patientID'] == patient_id].X))
    elif correlation_type == 'spearman':
        subset_data = adata[adata.obs['patientID'] == patient_id].X.toarray()
        corr_matrix, _ = spearmanr(subset_data, axis=0, nan_policy='omit')
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")
    
    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix)

    del corr_matrix
    gc.collect()

    process = psutil.Process(os.getpid())
    print(f"Memory usage inside patient {patient_id}: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bn.nanmean(rank_normalized_matrix)

    tracemalloc.stop_snapshot(f"After processing patient {patient_id}")
    tracemalloc.display_top()

    print(f"Finished processing patient {patient_id}")
    return patient_id, {
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }

def compute_gene_correlation(adata, correlation_type='pearson', replace_nans=True, 
                             min_cells_threshold=20, category=None, file_path="", 
                             output_dir=""):
    print("Starting correlation computation...")
    file_path = Path(file_path)
    output_dir = Path(output_dir)

    print(f"Processing subset for category: {category}")

    subset = adata[adata.obs['disease'] == category]
    patient_ids = subset.obs['patientID'].unique()
    total_patients = len(patient_ids)
    print(f"Total number of patients in {category}: {total_patients}")

    process = psutil.Process(os.getpid())
    print(f"Memory usage before loop: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    tracemalloc.start_snapshot("Before processing all patients")

    aggregated_matrix = None
    gene_names = None
    patient_count = 0

    for patient_id in patient_ids:
        print(f"Processing patient ID: {patient_id}")

        result = process_patient(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)

        if result is not None and result[1] is not None:
            if gene_names is None:
                gene_names = result[1]["gene_names"]

            if aggregated_matrix is None:
                aggregated_matrix = result[1]["rank_normalized_matrix"]
            else:
                aggregated_matrix += result[1]["rank_normalized_matrix"]

            patient_count += 1
            print(f"Patient {patient_id} matrix aggregated.")
            print(f"aggregated_matrix size after patient {patient_id}: {aggregated_matrix.shape}")
            print(f"Memory usage of aggregated_matrix after patient {patient_id}: {aggregated_matrix.nbytes / 1024 ** 2:.2f} MB")
            print(f"Memory usage after processing patient {patient_id}: {process.memory_info().rss / 1024 ** 2:.2f} MB")
        else:
            print(f"No valid matrix for patient {patient_id}.")
    
    if aggregated_matrix is not None and patient_count > 0:
        aggregated_matrix = aggregated_matrix / patient_count

    category_file_path = output_dir / f"{category}_aggregated_matrix.h5"
    with h5py.File(category_file_path, 'w') as f:
        f.create_dataset("aggregated_rank_normalized_matrix", data=aggregated_matrix.toarray())
        f.attrs["gene_names"] = np.array(gene_names, dtype="S")

    print(f"Results for {category} saved to {category_file_path}.")

    tracemalloc.stop_snapshot("After processing all patients")
    tracemalloc.display_top()

    print(f"Memory usage before cleanup: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    # Cleanup
    del subset
    del aggregated_matrix
    gc.collect()

    tracemalloc.start_snapshot("After cleanup")
    tracemalloc.display_top()

    print(f"Memory usage after cleanup: {process.memory_info().rss / 1024 ** 2:.2f} MB")
    print("Cleanup completed.")

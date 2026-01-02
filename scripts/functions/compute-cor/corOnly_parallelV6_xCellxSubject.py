
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
    Process data for a single patient: compute correlation and rank normalize.

    Returns:
        tuple -- Patient ID and a dictionary with processed edge lists and rank-normalized matrix.
    """
    #subset_data = adata[adata.obs['patientID'] == patient_id].X
    #subset_data = scipy.sparse.csr_matrix(subset_data)

    if adata[adata.obs['patientID'] == patient_id].X.shape[0] < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.")
        return patient_id, None
    
    # Compute correlation matrix
    if correlation_type == 'pearson':
        #corr_matrix = sparse_corr(subset_data)
        corr_matrix = sparse_corr(scipy.sparse.csr_matrix(adata[adata.obs['patientID'] == patient_id].X))
    elif correlation_type == 'spearman':
        corr_matrix, _ = spearmanr(subset_data, axis=0, nan_policy='omit')
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")
    
    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix) 

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bottleneck.nanmean(rank_normalized_matrix)

    print(f"Finished processing patient {patient_id}")
    del corr_matrix 
    #del subset_data
    return patient_id, {
        'rank_normalized_matrix': rank_normalized_matrix,
        'gene_names': adata.var_names.tolist()
    }


def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             category=None, 
                             file_path="", 
                             output_dir="",
                             xCellxSubject=False):
    """
    Compute gene-gene correlation matrices.
    
    Parameters
    ----------
    adata : AnnData
        Single-cell data.
    correlation_type : str
        'pearson' or 'spearman'.
    replace_nans : bool
        Whether to replace NaNs.
    min_cells_threshold : int
        Minimum number of cells to include a patient.
    category : str or None
        Disease or condition to subset. If 'ALL' or None, use all samples.
    xCellxSubject : bool
        If True, compute correlation across all patients at once; 
        if False, compute per-patient and aggregate.
    """
    print("Starting correlation computation...")
    file_path = Path(file_path)
    output_dir = Path(output_dir)

    # Handle category logic
    if category in [None, "ALL", "None", "none", ""]:
        subset = adata
        print("⚙️ Using all cells (no disease filtering).")
        category_label = "ALL"
    else:
        print(f"⚙️ Subsetting for category: {category}")
        subset = adata[adata.obs['disease'] == category].copy()
        category_label = category

    if subset.shape[0] == 0:
        print(f"⚠️ No cells found for category '{category_label}'. Skipping.")
        return

    gene_names = subset.var_names
    aggregated_matrix = None

    if xCellxSubject:
        print(f"🔹 Computing correlation across all patients for {category_label}")

        # Filter patients by min_cells_threshold
        patient_counts = subset.obs['patientID'].value_counts()
        valid_patients = patient_counts[patient_counts > min_cells_threshold].index
        filtered_subset = subset[subset.obs['patientID'].isin(valid_patients)]

        print(f"  Total patients after filtering: {len(valid_patients)}")
        print(f"  Total cells after filtering: {filtered_subset.shape[0]}")

        if filtered_subset.shape[0] == 0:
            print("⚠️ No patients with sufficient cells. Skipping correlation.")
            aggregated_matrix = None
            patient_count = 0
        else:
            expr_matrix = filtered_subset.X

            if correlation_type == 'pearson':
                aggregated_matrix = sparse_corr(scipy.sparse.csr_matrix(expr_matrix))
            elif correlation_type == 'spearman':
                aggregated_matrix, _ = spearmanr(expr_matrix, axis=0, nan_policy='omit')
            else:
                raise ValueError(f"Unsupported correlation_type: {correlation_type}")

            # Fill diagonal and handle NaNs
            np.fill_diagonal(aggregated_matrix, 1)
            aggregated_matrix = rank(aggregated_matrix)

            if replace_nans:
                aggregated_matrix[np.isnan(aggregated_matrix)] = bottleneck.nanmean(aggregated_matrix)

            patient_count = len(valid_patients)

    else:
        print(f"🔹 Computing correlation per patient for {category_label}")
        patient_ids = subset.obs['patientID'].unique()
        total_patients = len(patient_ids)
        print(f"  Total number of patients: {total_patients}")

        aggregated_matrix = np.zeros((subset.shape[1], subset.shape[1]))
        patient_count = 0

        for patient_id in patient_ids:
            print(f"  Processing patient ID: {patient_id}")
            result = process_patient(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)
            if result is not None and result[1] is not None:
                patient_matrix = result[1]["rank_normalized_matrix"]
                if aggregated_matrix is None:
                    aggregated_matrix = patient_matrix
                else:
                    aggregated_matrix.__iadd__(patient_matrix)
                patient_count += 1
                print(f"  ✅ Patient {patient_id} matrix aggregated.")
                del patient_matrix
                gc.collect()
            else:
                print(f"  ⚠️ No valid matrix for patient {patient_id}.")

        if aggregated_matrix is not None and patient_count > 0:
            aggregated_matrix = rank(aggregated_matrix)

    # Convert aggregated matrix to DataFrame
    aggregated_matrixEdge = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

    # Summary
    summary_df = pd.DataFrame({
        "category": [category_label],
        "total_genes": [len(gene_names)],
        "total_single_cells": [subset.shape[0]],
        "total_patients": [patient_count]
    })

    category_data = {
        'aggregated_rank_normalized_matrix': aggregated_matrix,
        'edgeList': aggregated_matrixEdge, 
        'gene_names': gene_names,
        'data_summary': summary_df, 
        'file_path': str(file_path).replace(".h5ad", f"_{category_label}.h5ad")
    }
    category_file_path = Path(category_data["file_path"])

    # Log memory usage
    process = psutil.Process(os.getpid())
    print(f"🧠 Memory usage before saving: {process.memory_info().rss / 1024 ** 2:.2f} MB")

    print(f"💾 Saving results for category: {category_label}")
    save_to_hdf5(category_data, output_dir, category_file_path)
    print(f"✅ Results for {category_label} saved successfully to {category_file_path}.")

    # Clean up memory
    del subset, aggregated_matrix, aggregated_matrixEdge, category_data, summary_df
    gc.collect()
    print(f"🧹 Memory usage after cleanup: {process.memory_info().rss / 1024 ** 2:.2f} MB")


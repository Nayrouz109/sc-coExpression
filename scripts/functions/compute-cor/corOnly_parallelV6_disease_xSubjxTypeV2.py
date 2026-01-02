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
import anndata as ad
from scipy import sparse


sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from saveCorDict import save_to_hdf5

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/DGE')
from create_PB import create_pseudo_bulk

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


def compute_gene_correlation(
    file_list, 
    correlation_type='pearson', 
    replace_nans=True, 
    category=None,     # e.g., ["CTL", "AD"] or "ALL" or None
    min_cells_threshold=20, 
    file_path="", 
    output_dir="",
    xType=True, 
    PB=True
):
    """
    Compute gene-gene correlation matrices, optionally stratified by disease category.
    Handles both 'ALL' and None gracefully by using all available samples.
    """

    def cpm_normalize(expr_df):
        """CPM normalize a pseudobulk expression matrix."""
        counts_sum = expr_df.sum(axis=1)
        return expr_df.div(counts_sum, axis=0) * 1e6

    # 🩵 Handle "ALL" or None categories
    if category is None or (isinstance(category, str) and category.upper() == "ALL"):
        category_label = "ALL"
        category_filter = None   # use all samples
    else:
        # Ensure list-like for consistency downstream
        category_label = category if isinstance(category, str) else "_".join(category)
        category_filter = [category] if isinstance(category, str) else category

    # ==========================================================
    # CASE 1: xType=True (aggregate across all files)
    # ==========================================================
    if PB and xType:
        pseudobulk_expr_all = []
        pseudobulk_meta_all = []

        for file in file_list:
            print(f"Processing {file} ...")
            adata = ad.read_h5ad(file)

            # 🔹 Subset to category if specified
            if category_filter is not None:
                subset = adata[adata.obs['disease'].isin(category_filter)].copy()
            else:
                subset = adata.copy()

            # 🔹 Filter patients with enough cells
            patient_counts = subset.obs['patientID'].value_counts()
            valid_patients = patient_counts[patient_counts > min_cells_threshold].index
            filtered_subset = subset[subset.obs['patientID'].isin(valid_patients)].copy()

            if filtered_subset.n_obs == 0:
                print(f"Skipping {file}: no valid patients after filtering.")
                continue

            # 🔹 Create pseudobulk and normalize
            pseudobulk_expr, pseudobulk_meta = create_pseudo_bulk(filtered_subset)
            pseudobulk_expr = cpm_normalize(pseudobulk_expr)

            pseudobulk_expr_all.append(pseudobulk_expr)
            pseudobulk_meta_all.append(pseudobulk_meta)

        # Concatenate across files
        if not pseudobulk_expr_all:
            print("No valid pseudobulk matrices found. Skipping correlation computation.")
            return

        expr_matrix = pd.concat(pseudobulk_expr_all, axis=0)
        gene_names = expr_matrix.columns

        # 🔹 Compute correlation
        if correlation_type == 'pearson':
            aggregated_matrix = sparse_corr(scipy.sparse.csr_matrix(expr_matrix))
        elif correlation_type == 'spearman':
            aggregated_matrix, _ = spearmanr(expr_matrix, axis=0, nan_policy='omit')
        else:
            raise ValueError(f"Unsupported correlation_type: {correlation_type}")
            
        # 🔹 Post-process correlation matrix
        np.fill_diagonal(aggregated_matrix, 1)
        aggregated_matrix = rank(aggregated_matrix)
        if replace_nans:
            aggregated_matrix[np.isnan(aggregated_matrix)] = bottleneck.nanmean(aggregated_matrix)

        aggregated_matrixEdge = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

        # 🔹 Prepare output
        category_file_path = Path(str(file_path).replace(".h5ad", f"_{category_label}.h5ad"))
        category_data = {
            'aggregated_rank_normalized_matrix': aggregated_matrix,
            'edgeList': aggregated_matrixEdge, 
            'gene_names': gene_names,
            'file_path': str(category_file_path)
        }

        print(f"Saving results for category: {category_label}")
        save_to_hdf5(category_data, output_dir, category_file_path)
        print(f"Results for {category_label} saved successfully to {category_file_path}.")

    # ==========================================================
    # CASE 2: xType=False — process each file independently
    # ==========================================================
    else:
        for file in file_list:
            print(f"\nProcessing {file} with xType=False ...")
            adata = ad.read_h5ad(file)

            # 🔹 Subset to category if specified
            if category_filter is not None:
                subset = adata[adata.obs['disease'].isin(category_filter)].copy()
            else:
                subset = adata.copy()

            # 🔹 Filter patients with enough cells
            patient_counts = subset.obs['patientID'].value_counts()
            valid_patients = patient_counts[patient_counts > min_cells_threshold].index
            filtered_subset = subset[subset.obs['patientID'].isin(valid_patients)].copy()

            if filtered_subset.n_obs == 0:
                print(f"Skipping {file}: no valid patients after filtering.")
                continue

            gene_names = filtered_subset.var_names

            # 🔹 Create pseudobulk or use raw expression
            if PB:
                pseudobulk_expr, pseudobulk_meta = create_pseudo_bulk(filtered_subset)
                pseudobulk_expr = cpm_normalize(pseudobulk_expr)
            else:
                pseudobulk_expr = pd.DataFrame(
                    filtered_subset.X.toarray() if hasattr(filtered_subset.X, "toarray") else filtered_subset.X,
                    index=filtered_subset.obs_names,
                    columns=gene_names
                )

            # 🔹 Compute correlation
            if correlation_type == 'pearson':
                aggregated_matrix = sparse_corr(scipy.sparse.csr_matrix(pseudobulk_expr))
            elif correlation_type == 'spearman':
                aggregated_matrix, _ = spearmanr(pseudobulk_expr, axis=0, nan_policy='omit')
            else:
                raise ValueError(f"Unsupported correlation_type: {correlation_type}")

            # 🔹 Post-process correlation matrix
            np.fill_diagonal(aggregated_matrix, 1)
            aggregated_matrix = rank(aggregated_matrix)
            if replace_nans:
                aggregated_matrix[np.isnan(aggregated_matrix)] = bottleneck.nanmean(aggregated_matrix)

            aggregated_matrixEdge = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

            # 🔹 Prepare output path
            save_file_path = Path(file).with_name(
                Path(file).stem + f"_{category_label}.h5ad"
            )

            category_data = {
                'aggregated_rank_normalized_matrix': aggregated_matrix,
                'edgeList': aggregated_matrixEdge,
                'gene_names': gene_names,
                'file_path': str(save_file_path)
            }

            print(f"Saving results for {file} — category: {category_label} ...")
            save_to_hdf5(category_data, output_dir, save_file_path)
            print(f"Results saved successfully to {save_file_path}.")

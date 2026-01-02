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
    category=None,     # list or vector of conditions e.g. ["CTL", "AD"]
    min_cells_threshold=20, 
    file_path = "", 
    output_dir="",
    xType=True, 
    PB=True
):


    def cpm_normalize(expr_df):
        """
        CPM normalize a pseudobulk expression matrix.
        """
        counts_sum = expr_df.sum(axis=1)
        return expr_df.div(counts_sum, axis=0) * 1e6

    pseudobulk_expr_all = []
    pseudobulk_meta_all = []

    if PB and xType:
        for file in file_list:
            print(f"Processing {file} ...")
            adata = ad.read_h5ad(file)

            # Ensure category is list-like
            if isinstance(category, str):category = [category]

            # Subset to condition(s)
            subset = adata[adata.obs['disease'].isin(category)].copy()

            # Filter patients with enough cells
            patient_counts = subset.obs['patientID'].value_counts()
            valid_patients = patient_counts[patient_counts > min_cells_threshold].index
            filtered_subset = subset[subset.obs['patientID'].isin(valid_patients)].copy()
            gene_names = filtered_subset.var_names
            aggregated_matrix = None

            # Skip if no valid patients
            if filtered_subset.n_obs == 0:
                print(f"Skipping {file}: no valid patients after filtering")
                continue

            # Create pseudobulk
            pseudobulk_expr, pseudobulk_meta = create_pseudo_bulk(filtered_subset)

            # CPM normalize
            pseudobulk_expr = cpm_normalize(pseudobulk_expr)

            # Collect
            pseudobulk_expr_all.append(pseudobulk_expr)
            pseudobulk_meta_all.append(pseudobulk_meta)

        # Concatenate across cell types
        pseudobulk_expr_concat = pd.concat(pseudobulk_expr_all, axis=0)
        #pseudobulk_meta_concat = pd.concat(pseudobulk_meta_all, axis=0)


        expr_matrix = pseudobulk_expr_concat

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

        # Convert aggregated matrix to DataFrame
        aggregated_matrixEdge = pd.DataFrame(aggregated_matrix, index=gene_names, columns=gene_names)

        category_data = {
        'aggregated_rank_normalized_matrix': aggregated_matrix,
        'edgeList': aggregated_matrixEdge, 
        'gene_names': gene_names,
        #'data_summary': summary_df, 
        'file_path': str(file_path).replace(".h5ad", f"_{category}.h5ad")
        }
        category_file_path = Path(category_data["file_path"])

        print(f"Saving results for category: {category}")
        save_to_hdf5(category_data, output_dir, category_file_path)
        print(f"Results for {category} saved successfully to {category_file_path}.")

    
    else:
        raise NotImplementedError("Currently only supports PB=True and xType=True case.")


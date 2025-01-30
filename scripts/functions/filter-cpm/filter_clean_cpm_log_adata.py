

import scipy.sparse
import anndata as ad
import numpy as np
import pandas as pd

def filter_and_clean_adata(data, selected_genes, sample_thr=0.02, cpm_normalize=False, log_transform=False):
    """
    1. Filters an AnnData object to include only selected genes, adds missing genes with zero expression, 
    2. filters single cells based on the number of expressed genes and a percentile threshold, 
    3. and optionally applies CPM normalization and log transformation.
    
    Parameters:
        data (AnnData): Input AnnData object.
        selected_genes (list): List of gene names to retain or add.
        sample_thr (float): Percentile threshold for filtering single cells based on the number of expressed genes.
        cpm_normalize (bool): Whether to apply CPM normalization to the filtered AnnData.
        log_transform (bool): Whether to apply log transformation to the CPM-normalized AnnData.

    Returns:
        AnnData: Processed AnnData object with selected genes and valid cells.
    """
    # Step 1: Filter single cells based on the number of expressed genes
    n_genes_per_cell = (data.X > 0).sum(axis=1).A1  # Count number of non-zero genes per cell
    perc_ranks = pd.Series(n_genes_per_cell).rank(pct=True)  # Percent rank for each cell
    valid_cells = perc_ranks > sample_thr  # Boolean mask for valid cells
    
    # Subset AnnData to include only valid cells
    data = data[valid_cells.values, :]

    # Step 2: Filter the AnnData for the selected genes present in the data
    current_genes = data.var_names.tolist()
    genes_to_keep = [gene for gene in selected_genes if gene in current_genes]
    data_filtered = data[:, genes_to_keep]
    data_filtered.X = scipy.sparse.csr_matrix(data_filtered.X)  # Ensure sparse format

    # Step 3: Identify missing genes
    missing_genes = [gene for gene in selected_genes if gene not in current_genes]

    # Step 4: Add missing genes with zero expression if any
    if missing_genes:
        missing_expr = scipy.sparse.csr_matrix((data_filtered.shape[0], len(missing_genes)))
        combined_matrix = scipy.sparse.hstack([missing_expr, data_filtered.X]).tocsc()

        # Create a new AnnData object with the updated data
        data_filtered_new = ad.AnnData(X=combined_matrix, obs=data_filtered.obs.copy())
        data_filtered_new.var_names = missing_genes + data_filtered.var_names.tolist()
    else:
        # No missing genes; return the filtered data directly
        data_filtered_new = data_filtered

    # Ensure sparse format is subscriptable
    data_filtered_new.X = scipy.sparse.csr_matrix(data_filtered_new.X)

    # Step 5: CPM Normalization (optional)
    if cpm_normalize:
        X = data_filtered_new.X
        if isinstance(X, scipy.sparse.spmatrix):  # Sparse matrix handling
            total_counts = X.sum(axis=1).A1  # Sum of counts per cell
        else:  # Dense matrix handling
            total_counts = X.sum(axis=1)
        total_counts[total_counts == 0] = 1  # Prevent division by zero
        data_filtered_new.X = (X.multiply(1e6 / total_counts[:, None]) 
                               if scipy.sparse.issparse(X) else X * 1e6 / total_counts[:, None])

    # Step 6: Log Transformation (optional)
    if log_transform:
        X = data_filtered_new.X
        if scipy.sparse.issparse(X):  # Sparse matrix handling
            data_filtered_new.X = X.log1p()
        else:  # Dense matrix handling
            data_filtered_new.X = np.log1p(X)

    # Ensure sparse format is subscriptable
    data_filtered_new.X = scipy.sparse.csr_matrix(data_filtered_new.X)

    return data_filtered_new







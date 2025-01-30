

import scipy.sparse
import anndata as ad
import numpy as np
import pandas as pd

def filter_and_clean_adata(data, selected_genes, sample_thr=0.02):
    """
    Filters an AnnData object to include only selected genes, adds missing genes with zero expression, 
    and filters single cells based on the number of expressed genes and a percentile threshold.
    
    Parameters:
        data (AnnData): Input AnnData object.
        selected_genes (list): List of gene names to retain or add.
        sample_thr (float): Percentile threshold for filtering single cells based on the number of expressed genes.

    Returns:
        AnnData: Filtered AnnData object with selected genes and valid cells.
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

    return data_filtered_new


## this is how to run the functon 
#import sys
#sys.path.append('/home/nelazzabi/rewiring/scripts/functions')

#from filter_clean_cpm-log_adata import filter_and_clean_adata

#filtered_data = filter_and_clean_adata(
    #data, 
    #selected_genes, 
    #sample_thr=0.02, 
    #cpm_normalize=True, 
    #log_transform=False
#)
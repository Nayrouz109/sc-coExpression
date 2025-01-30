

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

def clean_ctmat(adata, gene_thr=0.02, sample_thr=0.02):
    """
    Cleans an AnnData object by filtering out unwanted genes and cells.

    Parameters:
        adata (AnnData): The input AnnData object.
        gene_thr (float): The minimum fraction of cells in which a gene must be expressed to be retained.
        sample_thr (float): The minimum percentile of genes expressed per cell for a cell to be retained.

    Returns:
        AnnData: A new AnnData object with filtered genes and cells.
    """
    # Calculate the number of cells per gene and filter genes
    n_cells_tot = adata.n_obs
    n_cells_thr = gene_thr * n_cells_tot
    gene_counts = (adata.X > 0).sum(axis=0).A1 if isinstance(adata.X, csr_matrix) else (adata.X > 0).sum(axis=0)
    genes_to_keep = gene_counts > n_cells_thr


    # Filter genes
    adata = adata[:, genes_to_keep]

    # Calculate the number of genes per cell and filter cells
    gene_per_cell = np.array((adata.X > 0).sum(axis=1)).flatten()
    #gene_per_cell = (adata.X > 0).sum(axis=1).A1 if isinstance(adata.X, csr_matrix) else (adata.X > 0).sum(axis=1)
    #if not isinstance(gene_per_cell, np.ndarray):
        #gene_per_cell = np.asarray(gene_per_cell)
    cell_percentiles = pd.Series(gene_per_cell).rank(pct=True)
    cells_to_keep = cell_percentiles > sample_thr

    # Filter cells
    adata = adata[cells_to_keep, :]

    return adata




#cleaned_adata = clean_ctmat(data, gene_thr=0.02, sample_thr=0.02)

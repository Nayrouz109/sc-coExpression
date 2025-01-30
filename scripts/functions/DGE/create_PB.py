



### this is function to create PB matrices from single cell gene expression matrices 

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc 
from anndata import read_h5ad
from pathlib import Path


def create_pseudo_bulk(adata):
    """
    Create pseudobulk data by summing gene expression for each patient.

    Parameters:
    - adata: AnnData object containing single-cell data
    - patient_col: Column name in adata.obs indicating patient IDs

    Returns:
    - pseudobulk_expr: DataFrame with pseudobulk expression data (genes x patients)
    - pseudobulk_meta: DataFrame with metadata for pseudobulk samples
    """
    expr = adata.X  # Gene expression matrix (sparse or dense)
    meta = adata.obs

    # Ensure expression is dense for aggregation
    if not isinstance(expr, np.ndarray):
        expr = expr.toarray()

    # Group cells by patient and sum the expression values
    meta['patientID'] = meta['patientID'].astype(str).fillna('unknown')
    patient_groups = meta.groupby("patientID")
    pseudobulk_expr = patient_groups.apply(
        lambda group: expr[group.index.to_series().map(meta.index.get_loc), :].sum(axis=0)
    )
    pseudobulk_expr.index.name = None
    pseudobulk_expr = pd.DataFrame(
        np.vstack(pseudobulk_expr.values),
        index=pseudobulk_expr.index,
        columns=adata.var.index
    )
    # Metadata for pseudobulk samples
    pseudobulk_meta = meta["patientID"].drop_duplicates().reset_index(drop=True)
    pseudobulk_meta.index = pseudobulk_expr.index

    return pseudobulk_expr, pseudobulk_meta





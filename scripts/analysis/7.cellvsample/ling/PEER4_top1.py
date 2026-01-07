#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

"""
Extract gene-gene correlation values for PEER4-selected (SNAP) genes

Also extracts:
- top1: top 1% strongest coexpression edges genome-wide (per matrix)
- background_random1000: 1000 random gene pairs genome-wide (per matrix)

Author: Nairuz
Date: 2025-10-08
"""

import os
import argparse
import pandas as pd
import anndata as ad
import numpy as np
from pathlib import Path


# -----------------------------------
# Extract top genes by cell type and loading sign
# -----------------------------------
def extract_gene_set(
    data: pd.DataFrame,
    top_n: int,
    cell_type=None,
    loading_sign: str = None
) -> pd.DataFrame:
    df = data.copy()

    # Filter by cell type if provided
    if cell_type is not None:
        if isinstance(cell_type, str):
            cell_type = [cell_type]
        df = df.query("cell_type in @cell_type")

    # Determine ranking mode
    if loading_sign == "positive":
        filtered = df.query("sign_PEER_4 == 1").sort_values("PEER_4", ascending=False)
        msg = f"positive genes for {cell_type}"
    elif loading_sign == "negative":
        filtered = df.query("sign_PEER_4 == -1").sort_values("PEER_4", ascending=True)
        msg = f"negative genes for {cell_type}"
    elif loading_sign is None:
        df["abs_PEER_4"] = df["PEER_4"].abs()
        filtered = df.sort_values("abs_PEER_4", ascending=False)
        msg = "genes by |PEER_4| (agnostic of cell type and sign)"
    else:
        raise ValueError("loading_sign must be 'positive', 'negative', or None")

    top_genes = filtered.head(top_n)
    print(f" Extracted top {len(top_genes)} {msg}\n")
    return top_genes


def extract_top1_edges(corr_values: np.ndarray, genes: pd.Index, cell_type: str, study: str) -> pd.DataFrame:
    """
    Extract top 1% strongest edges from a correlation matrix (upper triangle, k=1),
    ranked by absolute correlation magnitude.
    """
    n_genes = len(genes)
    if n_genes < 2:
        return pd.DataFrame()

    iu = np.triu_indices(n_genes, k=1)
    vals = corr_values[iu]  # correlations for unique pairs
    m = len(vals)
    if m == 0:
        return pd.DataFrame()

    k = max(1, int(np.floor(0.01 * m)))  # top 1%
    abs_vals = np.abs(vals)

    # indices of top-k absolute correlations (unordered)
    topk_idx = np.argpartition(abs_vals, -k)[-k:]

    g1 = genes[iu[0][topk_idx]]
    g2 = genes[iu[1][topk_idx]]
    cor = vals[topk_idx]

    return pd.DataFrame({
        "gene1": g1,
        "gene2": g2,
        "correlation": cor,
        "cell_type": cell_type,
        "gene_list": "top1",
        "study": study
    })


def main(input_dir, output_dir, suffix=""):
    os.makedirs(output_dir, exist_ok=True)

    # ---------------------------
    # Load PEER loadings and ribosomal genes
    # ---------------------------
    peer_file = "/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/lings_PEER4_gene_loadings_for_nairuz.csv"
    ribo_file = "/home/nelazzabi/rewiring/data/genes/ribosomal/group-1054.csv"

    print("\n Reading PEER loadings file...")
    genes2 = pd.read_csv(peer_file)

    print(" Reading ribosomal genes file...")
    ribo = pd.read_csv(ribo_file, skiprows=1)
    ribo = ribo[ribo["Group"].isin(["L ribosomal proteins", "S ribosomal proteins"])]

    # ---------------------------
    # Extract gene sets
    # ---------------------------
    SNAP = extract_gene_set(genes2, 1000, cell_type=None, loading_sign=None)

    gene_sets = {
        "SNAP": SNAP["gene"].tolist(),
        "ribo": ribo["Approved symbol"].tolist(),
    }

    # ---------------------------
    # Define dataset patterns
    # ---------------------------
    input_dir = Path(input_dir)
    patterns = {
        "ling-CTL": "*ling*CTL*.h5ad",
        "ROSMAP-Mathys-2019-ALL": "*ROSMAP-Mathys-2019*ALL*.h5ad",
        "Lau-2020-ALL": "*Lau-2020*ALL*.h5ad",
        "Lim-2022-ALL": "*Lim-2022*ALL*.h5ad",
        "Nagy-2020-ALL": "*Nagy-2020*ALL*.h5ad",
        "Pineda-2021-ALL": "*Pineda-2021*ALL*.h5ad",
        "Velmeshev-2019-ALL": "*Velmeshev-2019*ALL*.h5ad",
    }

    # ---------------------------
    # Iterate per pattern (study)
    # ---------------------------
    for study, pattern in patterns.items():
        study_files = list(input_dir.glob(pattern))
        if not study_files:
            print(f" No files found for pattern: {pattern}")
            continue

        print(f"\n Processing study: {study} ({len(study_files)} file(s) found)")
        results = []

        for file_path in study_files:
            ct = file_path.stem.split("_")[0]
            print(f" Processing {ct}...")

            adata = ad.read_h5ad(file_path)
            corr_matrix = adata.X
            genes = adata.var_names

            # Convert to dense numpy array

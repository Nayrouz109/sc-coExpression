#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

"""
Extract gene-gene correlation values for PEER4-selected genes
from Ling et al. (2024) single-cell coexpression datasets.

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
# Step 2: Extract top genes by cell type and loading sign
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

    # Determine which ranking mode to use
    if loading_sign == "positive":
        filtered = df.query("sign_PEER_4 == 1").sort_values("PEER_4", ascending=False)
        msg = f"positive genes for {cell_type}"
    elif loading_sign == "negative":
        filtered = df.query("sign_PEER_4 == -1").sort_values("PEER_4", ascending=True)
        msg = f"negative genes for {cell_type}"
    elif loading_sign is None:
        # Agnostic to sign and cell type: use absolute PEER_4
        df["abs_PEER_4"] = df["PEER_4"].abs()
        filtered = df.sort_values("abs_PEER_4", ascending=False)
        msg = "genes by |PEER_4| (agnostic of cell type and sign)"
    else:
        raise ValueError("loading_sign must be 'positive', 'negative', or None")

    top_genes = filtered.head(top_n)
    print(f"🧬 Extracted top {len(top_genes)} {msg}\n")
    return top_genes



# -----------------------------------
# Main Execution
# -----------------------------------
def main(input_dir, output_dir,  suffix=""):
    os.makedirs(output_dir, exist_ok=True)
    peer_file = "/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/lings_PEER4_gene_loadings_for_nairuz.csv"
    #peer_file = "/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/Ling/mine_PEER4_gene_loadings_for_nairuz.csv"
    ribo = pd.read_csv("/home/nelazzabi/rewiring/data/genes/ribosomal/group-1054.csv")
    ribo = pd.read_csv("/home/nelazzabi/rewiring/data/genes/ribosomal/group-1054.csv", skiprows=1)
    ribo = ribo[ribo["Group"].isin(["L ribosomal proteins", "S ribosomal proteins"])]

    print("\n📥 Reading PEER loadings file...")
    genes2 = pd.read_csv(peer_file)

    # Extract top genes for each group
    astro_pos = extract_gene_set(genes2, 1000, "astrocyte", "positive")
    astro_neg = extract_gene_set(genes2, 1000, "astrocyte", "negative")
    neur_pos  = extract_gene_set(genes2, 1000, ["glutamatergic", "gabaergic"], "positive")
    neur_neg  = extract_gene_set(genes2, 1000, ["glutamatergic", "gabaergic"], "negative")
    top_abs   = extract_gene_set(genes2, 1000, cell_type=None, loading_sign=None)

    gene_sets = {
        "astro_pos": astro_pos["gene"].tolist(),
        "astro_neg": astro_neg["gene"].tolist(),
        "neur_pos": neur_pos["gene"].tolist(),
        "neur_neg": neur_neg["gene"].tolist(),
        "top_abs": top_abs["gene"].tolist(),
        "ribo": ribo["Approved symbol"].tolist()  # add ribo gene set
    }

    # Find Ling-2024 files
    input_dir = Path(input_dir)
    pattern = "*ling*CTL*.h5ad"
    ling_files = list(input_dir.glob(pattern))
    results = []

    # Iterate through Ling files and extract correlations
    # Iterate through Ling-CTL files and extract correlations
    for file_path in ling_files:
        # Extract cell type from filename (before the first underscore)
        ct = file_path.stem.split("_")[0]
        print(f"🔍 Processing {ct}...")

        adata = ad.read_h5ad(file_path)
        corr_matrix = adata.X
        genes = adata.var_names

        # Convert to dense if sparse
        print(f" converting to a dataframe...")
        corr_df = pd.DataFrame(
            corr_matrix.todense() if hasattr(corr_matrix, "todense") else corr_matrix,
            index=genes,
            columns=genes
        )

        # Loop over gene sets
        print(f" looping over the genes...")
        for gene_list_name, gene_list in gene_sets.items():
            valid_genes = [g for g in gene_list if g in genes]
            if len(valid_genes) < 2:
                print(f"⚠️ Not enough valid genes from {gene_list_name} in {ct}")
                continue

            df = (
                corr_df.loc[valid_genes, valid_genes]
                .stack()
                .reset_index()
                .rename(columns={"level_0": "gene1", "level_1": "gene2", 0: "correlation"})
                .assign(cell_type=ct, gene_list=gene_list_name)
            )

            # Remove self-correlations and duplicate pairs
            df = df[df["gene1"] != df["gene2"]]
            df["pair"] = df.apply(lambda x: tuple(sorted([x["gene1"], x["gene2"]])), axis=1)
            df = df.drop_duplicates(subset=["cell_type", "pair"]).drop(columns="pair")

            results.append(df)

    # Combine all results
    combined_df = pd.concat(results, ignore_index=True)

    # Save to CSV
    #filename = f"ruiPEER4_gene_correlations_all_sets{suffix}.csv"
    filename = f"PEER4_gene_correlations_all_sets{suffix}.csv"
    output_csv = os.path.join(output_dir, filename)
    combined_df.to_csv(output_csv, index=False)

    print(f"\n✅ Saved combined correlation data (with gene_list and cell_type) to:\n{output_csv}")
    print("🎉 Finished successfully!\n")


# -----------------------------------
# CLI Entry Point
# -----------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract PEER4-related gene correlations for Ling-2024 coexpression matrices."
    )
    parser.add_argument("--input_dir", required=True, help="Path to input directory containing .h5ad files")
    parser.add_argument("--output_dir", required=True, help="Path to save results")
    parser.add_argument("--suffix", default="", help="Optional suffix to add to the output CSV filename")

    args = parser.parse_args()

    main(args.input_dir, args.output_dir, args.suffix)


#python extract_peer4_corrs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/agg-mtx \ 
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4 \ 
#  --suffix _xCellxSubject


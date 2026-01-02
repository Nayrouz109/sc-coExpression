#!/~/anaconda3/envs/homl3/bin/python
import os
import gzip
import pandas as pd
import scanpy as sc
from scipy.io import mmread

def read_mtx_gz(file_path):
    """Read a gzipped matrix.mtx file into a sparse matrix."""
    with gzip.open(file_path, 'rb') as f:
        return mmread(f).tocsr()

def create_anndata_per_celltype(input_dir, output_dir, celltype_file):
    # Mapping from raw cell type names to canonical names
    celltype_map = {
        "Excit": "Exc",
        "Inhit": "Inh",
        "Oligo": "Oli",
        "Astro": "Ast",
        "Mic": "Mic",
        "Opc": "Opc"
    }
    valid_cell_types = list(celltype_map.values())

    # Read cell type annotations
    celltype_df = pd.read_excel(celltype_file)
    celltype_df['ID'] = celltype_df['ID'].astype(str)

    # Map raw names → canonical names
    celltype_df['Cell_type'] = celltype_df['Cell_type'].map(celltype_map)

    # Keep only valid cell types
    celltype_df = celltype_df[celltype_df['Cell_type'].notna()]

    # Collect all matrix files
    mtx_files = [f for f in os.listdir(input_dir) if f.endswith("_matrix.mtx.gz")]

    for mtx_file in mtx_files:
        prefix = mtx_file.replace("_matrix.mtx.gz", "")
        barcodes_file = os.path.join(input_dir, prefix + "_barcodes.tsv.gz")
        features_file = os.path.join(input_dir, prefix + "_features.tsv.gz")
        mtx_path = os.path.join(input_dir, mtx_file)

        # Read files
        print(f"Processing {prefix} ...")
        X = read_mtx_gz(mtx_path)
        barcodes = pd.read_csv(barcodes_file, header=None, names=['barcode'])
        features = pd.read_csv(features_file, sep="\t", header=None)

        # Fix barcodes suffix: -1 → _1
        barcodes['barcode'] = barcodes['barcode'].str.replace(r'-1$', '_1', regex=True)

        # Make sure matrix shape matches
        assert X.shape[0] == features.shape[0], f"{prefix}: feature length mismatch"
        assert X.shape[1] == barcodes.shape[0], f"{prefix}: barcode length mismatch"

        # Create AnnData (cells x genes)
        adata = sc.AnnData(X.T)
        adata.var['gene_ids'] = features.iloc[:, 0].values
        adata.var['gene_names'] = features.iloc[:, 1].values
        adata.obs['cell_ids'] = barcodes['barcode'].values
        adata.obs_names = barcodes['barcode'].values
        adata.var_names = features.iloc[:, 1].values

        # Merge with cell type annotations
        # Make sure ID is the index for lookup
        celltype_df = celltype_df.set_index("ID")

        # Map cell type onto obs using the index
        adata.obs["Cell_type"] = adata.obs.index.map(celltype_df["Cell_type"])

        # Split by canonical cell types and save
        for cell_type in valid_cell_types:
            adata_sub = adata[adata.obs['Cell_type'] == cell_type].copy()
            if adata_sub.n_obs == 0:
                continue
            outfile = os.path.join(output_dir, f"{cell_type}_Lau-2020.h5ad")
            adata_sub.write(outfile)
            print(f"Saved {outfile} ({adata_sub.n_obs} cells, {adata_sub.n_vars} genes)")

if __name__ == "__main__":
    input_dir = "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/lau"
    output_dir = "/cosmos/data/project-data/NW-rewiring/data/pans/0.prcsd-raw"
    celltype_file = "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/lau/Lau_et_al_PNAS_AD_snRNAseq_celltype_subcluster_ID.xlsx"

    os.makedirs(output_dir, exist_ok=True)
    create_anndata_per_celltype(input_dir, output_dir, celltype_file)

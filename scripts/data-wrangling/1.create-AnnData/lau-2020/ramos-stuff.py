

import scanpy as sc
import anndata as ad
import os
import re
import pandas as pd

# Path to your raw directory
raw_dir = "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/ramos/GSE217511_RAW"

# Collect triplets of files
matrix_files = sorted([f for f in os.listdir(raw_dir) if f.endswith("matrix.mtx.gz")])

all_adatas = []

for mtx_file in matrix_files:
    # Derive sample prefix (e.g., GSM6720852_9C)
    prefix = mtx_file.replace("_matrix.mtx.gz", "")
    barcodes_file = os.path.join(raw_dir, prefix + "_barcodes.tsv.gz")
    features_file = os.path.join(raw_dir, prefix + "_features.tsv.gz")
    mtx_file = os.path.join(raw_dir, mtx_file)

    # Extract sample/region ID from prefix (after GSM number -> e.g. "9C")
    match = re.search(r"_(\d+[A-Z])$", prefix)
    sample_id = match.group(1) if match else prefix

    # --- Read files explicitly ---
    X = sc.read_mtx(mtx_file).X.tocsr()  # expression matrix (sparse)
    barcodes = pd.read_csv(barcodes_file, header=None)[0].astype(str)
    features = pd.read_csv(features_file, header=None, sep="\t")

    # Format barcodes with sample suffix
    barcodes = barcodes + "_" + sample_id

    # Build var dataframe (gene-level metadata)
    var = pd.DataFrame({
        "gene_id": features[0].values,
        "gene_name": features[1].values,
        "feature_type": features[2].values
    }, index=features[1].values)   # index = gene symbols

    # --- Build obs dataframe (cell-level metadata) ---
    obs = pd.DataFrame(index=barcodes)
    obs["sample_id"] = sample_id

    # --- Create AnnData object ---
    adata = ad.AnnData(X=X.T, obs=obs, var=var)  
    adata.var_names_make_unique()

    all_adatas.append(adata)

# Concatenate across all samples
combined = ad.concat(all_adatas, join="outer", label="sample_id", fill_value=0)
combined.obs = combined.obs.reset_index().rename(columns={"index": "cell_id"}) #assign a cell ID column
combined.obs = combined.obs.rename(columns={combined.obs.columns[0]: "cell_id"})

combined.obs['sample'] = combined.obs['cell_id'].str.extract(r'_([A-Za-z0-9]+)$')

metadata_file = "/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/ramos/clean_metadata.csv"
metadata = pd.read_csv(metadata_file)
metadata = metadata.rename(columns={"cellids": "cell_id"})



# Number of cells in meta that are in combined.obs
metadata["cellids"].isin(combined.obs["cell_id"]).sum()








combined.obs = combined.obs.merge(metadata[["cell_id", "celltypes"]],
                                  on="cell_id",
                                  how="left")



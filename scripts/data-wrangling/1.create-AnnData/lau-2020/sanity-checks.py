




import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import anndata as ad 
import os

# Set working directory
data_dir = "/cosmos/data/project-data/NW-rewiring/data/pans/0.prcsd-raw"

# List of .h5ad files
files = [
    "Lau-2020.h5ad",
    "Lim-2022.h5ad",
    "Nagy-2020.h5ad",
    "Pineda-2021.h5ad",
    "Velmeshev-2019.h5ad",
]

for fname in files:
    fpath = os.path.join(data_dir, fname)
    print(f"\n=== {fname} ===")
    
    # Load in backed mode first
    adata = sc.read_h5ad(fpath, backed="r")
    
    # Check if per-cell QC already exists
    if "n_counts" in adata.obs.columns:
        counts = adata.obs["n_counts"].values
    else:
        print("⚠️  n_counts not found, loading full data (may be slow for big files)...")
        adata_full = sc.read_h5ad(fpath)  # loads into memory
        counts = np.array(adata_full.X.sum(axis=1)).flatten()
    
    print(f"Cells: {counts.shape[0]}")
    print(f"Mean counts per cell: {np.mean(counts):.2f}")
    print(f"Median counts per cell: {np.median(counts):.2f}")
    print(f"Min: {np.min(counts)}, Max: {np.max(counts)}")
    
    # Histogram
    plt.figure(figsize=(5,3))
    plt.hist(counts, bins=50, color="steelblue")
    plt.title(f"{fname} - counts per cell")
    plt.xlabel("Total counts per cell")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.tight_layout()
    plt.show()






#########################################################
#### summary table 
#########################################################

import os
import scanpy as sc
import numpy as np
import pandas as pd

data_dir = "/cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/raw"
files = [f for f in os.listdir(data_dir) if f.endswith(".h5ad")]

summary = []

for f in sorted(files):
    file_path = os.path.join(data_dir, f)
    adata = sc.read_h5ad(file_path)
    
    # Compute total counts per cell
    counts = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
    
    summary.append({
        "File": f,
        "Cells": counts.shape[0],
        "Mean_counts_per_cell": np.mean(counts),
        "Median_counts_per_cell": np.median(counts),
        "Min_counts": np.min(counts),
        "Max_counts": np.max(counts)
    })

# Convert to DataFrame
summary_df = pd.DataFrame(summary)

# Display nicely
print(summary_df)

#########################################################
#########################################################








og = ad.read_h5ad("/cosmos/data/project-data/NW-rewiring/data/pans/0.prcsd-raw/Nagy-2020.h5ad")
adata = og 
counts = np.array(adata.X.sum(axis=1)).flatten()
print(f"Cells: {counts.shape[0]}")
print(f"Mean counts per cell: {np.mean(counts):.2f}")
print(f"Median counts per cell: {np.median(counts):.2f}")
print(f"Min: {np.min(counts)}, Max: {np.max(counts)}")


og1 = ad.read_h5ad("/cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/cpm/Mic_Nagy-2020.h5ad")
counts = np.array(og1.X.sum(axis=1)).flatten()
print(f"Cells: {counts.shape[0]}")
print(f"Mean counts per cell: {np.mean(counts):.2f}")
print(f"Median counts per cell: {np.median(counts):.2f}")
print(f"Min: {np.min(counts)}, Max: {np.max(counts)}")
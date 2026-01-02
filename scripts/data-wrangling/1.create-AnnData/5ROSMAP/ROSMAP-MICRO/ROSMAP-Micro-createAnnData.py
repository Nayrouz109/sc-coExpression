
# this script was adopted from ling-2024 

import pandas as pd
import scipy.io
import anndata as ad
import os
from scipy.sparse import csr_matrix

path  = "/space/scratch/bxu_project/ROSMAP-Micro/" 
output_dir = "/space/scratch/bxu_project/ROSMAP-Micro/"
cell_type = "Mic"

count = "ROSMAP-Micro_counts.mtx"
cell_id = "ROSMAP-Micro_cell_ids.csv"
gene = "ROSMAP-Micro_genes.csv"
meta = "ROSMAP-Micro_meta_data.csv"




    
print(f"  Loading counts, cells, and metadata for {cell_type}...")
        
counts = scipy.io.mmread(os.path.join(path, count)).tocsc()
print(f"  Loaded counts matrix for {cell_type}")

cells = pd.read_csv(os.path.join(path, cell_id))
obs = pd.read_csv(os.path.join(path, meta)) 
genes = pd.read_csv(os.path.join(path, gene))
counts = counts.transpose()  # make it cells x genes
print(f"  Loaded cell IDs, metadata, and genes for {cell_type}")

    # Ensure row and column names match
    #counts = counts.transpose()  # make it cells x genes
genes.columns = ["gene_names"]
cells.columns = ["cell_ids"]
print(f"  Prepared data shapes: counts ({counts.shape}), cells ({cells.shape}), genes ({genes.shape})")

    # Create AnnData object
adata = ad.AnnData(X=counts, obs=obs, var=genes)
adata.obs_names = cells["cell_ids"]
adata.var_names = genes["gene_names"]
print(f"  Created AnnData object for {cell_type} with shape {adata.shape}")

    # Save to output directory
output_path = os.path.join(output_dir, f"{cell_type}_ROSMAP-Micro.h5ad")
adata.write_h5ad(output_path, compression="gzip")
print(f"  Successfully saved AnnData for {cell_type} to {output_path}")

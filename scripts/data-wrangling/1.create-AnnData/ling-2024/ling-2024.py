#!/~/anaconda3/envs/homl3/bin/python

### 00README 
# ling-2024 was previously processed into cell type rds 
# this script combines the extracted rds components into an AnnData format 

import pandas as pd
import scipy.io
import anndata as ad
import os
from scipy.sparse import csr_matrix

# Paths to rds components 
path  = "/cosmos/data/downloaded-data/AD-meta/ling-2024/data/extracted-rds-comps"
output_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw/"

# Define input files and output directory
file_info = {
    #"Ast": ["Ast_ling-2024_countsV2.mtx", "Ast_ling-2024_cell_ids.csv", "Ast_ling-2024_genes.csv", "Ast_ling-2024_meta_data.csv"],
    "Exc": [["Exc_batch1_ling-2024_countsV2.mtx", "Exc_batch2_ling-2024_countsV2.mtx"],
            ["Exc_batch1_ling-2024_cell_ids.csv", "Exc_batch2_ling-2024_cell_ids.csv"],
            ["Exc_batch1_ling-2024_genes.csv"], 
            ["Exc_batch1_ling-2024_meta_data.csv", "Exc_batch2_ling-2024_meta_data.csv"]]#,
    #"Inh": ["Inh_ling-2024_countsV2.mtx", "Inh_ling-2024_cell_ids.csv", "Inh_ling-2024_genes.csv", "Inh_ling-2024_meta_data.csv"],
    #"Mic": ["Mic_ling-2024_countsV2.mtx", "Mic_ling-2024_cell_ids.csv", "Mic_ling-2024_genes.csv", "Mic_ling-2024_meta_data.csv"],
    #"Oli": ["Oli_ling-2024_countsV2.mtx", "Oli_ling-2024_cell_ids.csv", "Oli_ling-2024_genes.csv", "Oli_ling-2024_meta_data.csv"],
    #"Opc": ["Opc_ling-2024_countsV2.mtx", "Opc_ling-2024_cell_ids.csv", "Opc_ling-2024_genes.csv", "Opc_ling-2024_meta_data.csv"]
}

def create_anndata(cell_type, files):
    print(f"\nProcessing cell type: {cell_type}")
    
    if cell_type == "Exc":  # Special handling for "Exc" to combine batch1 and batch2
        print(f"  Loading and concatenating counts, cells, and metadata for {cell_type}...")
        
        counts = [scipy.io.mmread(os.path.join(path, f)) for f in files[0]]
        #counts = scipy.sparse.vstack(counts)  # Concatenate along rows
        transposed_matrices = [matrix.transpose() for matrix in counts]
        counts = scipy.sparse.vstack(transposed_matrices).tocsc() 
        print(f"  Loaded and concatenated counts matrices for {cell_type}")

        cells = pd.concat([pd.read_csv(os.path.join(path, f)) for f in files[1]], ignore_index=True)
        obs = pd.concat([pd.read_csv(os.path.join(path, f)) for f in files[3]], ignore_index=True)
        print(f"  Loaded and concatenated cell IDs and metadata for {cell_type}")

        genes = pd.read_csv(os.path.join(path, files[2][0]))
        print(f"  Loaded genes file for {cell_type}")

    else:
        print(f"  Loading counts, cells, and metadata for {cell_type}...")
        
        counts = scipy.io.mmread(os.path.join(path, files[0])).tocsc()
        print(f"  Loaded counts matrix for {cell_type}")

        cells = pd.read_csv(os.path.join(path, files[1]))
        obs = pd.read_csv(os.path.join(path, files[3])) 
        genes = pd.read_csv(os.path.join(path, files[2]))
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
    output_path = os.path.join(output_dir, f"{cell_type}_ling-2024.h5ad")
    adata.write_h5ad(output_path, compression="gzip")
    print(f"  Successfully saved AnnData for {cell_type} to {output_path}")

# Create and save AnnData for each cell type
for cell_type, files in file_info.items():
    create_anndata(cell_type, files)

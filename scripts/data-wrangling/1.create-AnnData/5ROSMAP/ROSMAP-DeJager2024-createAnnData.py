
#!/~/anaconda3/envs/homl3/bin/python

# this script was adopted from ling-2024 

import pandas as pd
import scipy.io
import anndata as ad
import os
from scipy.sparse import csr_matrix

path  = "/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager/data/0-source/expr/rds/rds-comps/" 
output_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw"


file_info = {
    "Ast": ["Ast_ROSMAP-DeJager-2024_counts.mtx", "Ast_ROSMAP-DeJager-2024_cell_ids.csv", "Ast_ROSMAP-DeJager-2024_genes.csv", "Ast_ROSMAP-DeJager-2024_meta_data.csv"],
    "Exc": [["Exc1_ROSMAP-DeJager-2024_counts.mtx", "Exc2_ROSMAP-DeJager-2024_counts.mtx"],
            ["Exc1_ROSMAP-DeJager-2024_cell_ids.csv", "Exc2_ROSMAP-DeJager-2024_cell_ids.csv"],
            ["Exc1_ROSMAP-DeJager-2024_genes.csv"], 
            ["Exc1_ROSMAP-DeJager-2024_meta_data.csv", "Exc2_ROSMAP-DeJager-2024_meta_data.csv"]],
    "Inh": ["Inh_ROSMAP-DeJager-2024_counts.mtx", "Inh_ROSMAP-DeJager-2024_cell_ids.csv", "Inh_ROSMAP-DeJager-2024_genes.csv", "Inh_ROSMAP-DeJager-2024_meta_data.csv"],
    "Mic": ["Mic_ROSMAP-DeJager-2024_counts.mtx", "Mic_ROSMAP-DeJager-2024_cell_ids.csv", "Mic_ROSMAP-DeJager-2024_genes.csv", "Mic_ROSMAP-DeJager-2024_meta_data.csv"],
    "Oli": ["Oli_ROSMAP-DeJager-2024_counts.mtx", "Oli_ROSMAP-DeJager-2024_cell_ids.csv", "Oli_ROSMAP-DeJager-2024_genes.csv", "Oli_ROSMAP-DeJager-2024_meta_data.csv"]
}



def create_anndata(cell_type, files):
    print(f"\nProcessing cell type: {cell_type}")
    
    if cell_type == "Exc": 
        print(f"  Loading and concatenating counts, cells, and metadata for {cell_type}...")
        
        counts = [scipy.io.mmread(os.path.join(path, f)) for f in files[0]]
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
    output_path = os.path.join(output_dir, f"{cell_type}_ROSMAP-DeJager-2024.h5ad")
    adata.write_h5ad(output_path, compression="gzip")
    print(f"  Successfully saved AnnData for {cell_type} to {output_path}")

# Create and save AnnData for each cell type
for cell_type, files in file_info.items():
    create_anndata(cell_type, files)




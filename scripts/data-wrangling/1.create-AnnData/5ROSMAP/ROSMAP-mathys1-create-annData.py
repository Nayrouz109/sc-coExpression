#!/~/anaconda3/envs/homl3/bin/python


import pandas as pd
import scipy.io
import anndata as ad
import os
from scipy.sparse import csr_matrix


path  = "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/expr/rds-comps/"
output_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw"

all_files = sorted(os.listdir(path))
valid_prefixes = {"Ast", "Exc", "Mic", "Opc", "Oli", "Inh"}


for prefix in set(f.split("_")[0] + "_" + f.split("_")[1] for f in all_files):
    if prefix.split("_")[0] in valid_prefixes:
        # Define file paths based on the naming conventions
        count_file = os.path.join(path, f"{prefix}_counts.mtx")
        cell_file = os.path.join(path, f"{prefix}_cell_ids.csv")
        gene_file = os.path.join(path, f"{prefix}_genes.csv")
        meta_file = os.path.join(path, f"{prefix}_meta_data.csv")
    
        # Check if all files exist
        if all(os.path.exists(f) for f in [count_file, cell_file, gene_file, meta_file]):
            # Read in the data
            print(f"reading data for {prefix}")
            counts = scipy.io.mmread(count_file).tocsc()  # Sparse matrix format
            cells = pd.read_csv(cell_file)
            genes = pd.read_csv(gene_file)
            meta_data = pd.read_csv(meta_file)

            counts = counts.transpose()  # make it cells x genes
            genes.columns = ["gene_names"]
            cells.columns = ["cell_ids"]
        
            # Ensure dimensions match for AnnData creation 
            if counts.shape[1] == genes.shape[0] and counts.shape[0] == cells.shape[0]:
                # Create AnnData object
                print(f"creating annData for {prefix}")
                adata = ad.AnnData(X=counts, obs=meta_data, var=genes)
                adata.obs_names = cells["cell_ids"]
                adata.var_names = genes["gene_names"]
            
                # Store the AnnData object
                output_filename = os.path.join(output_dir, f"{prefix}.h5ad")
                adata.write_h5ad(output_filename, compression="gzip")

                print(f"saved AnnData for {prefix}")
            else:
                print(f"Dimension mismatch for {prefix}")
        else:
            print(f"Files missing for {prefix}")
    else:
        print(f"Skipping {prefix} as it is not in the valid prefixes list.")




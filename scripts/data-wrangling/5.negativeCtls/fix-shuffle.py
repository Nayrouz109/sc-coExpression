#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import os
import anndata as ad
import pandas as pd

# Directory containing the .h5ad files
data_dir = "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/01shuffledDisease/AD-Ctl"

# Loop over all .h5ad files in the directory
for filename in os.listdir(data_dir):
    if filename.endswith(".h5ad"):  # Ensure it's an .h5ad file
        file_path = os.path.join(data_dir, filename)
        
        # Load the .h5ad file
        adata = ad.read_h5ad(file_path)
        
        # Modify column names
        adata.obs = adata.obs.rename(columns={
            "disease": "old_disease",
            "random_group": "disease"
        })
        
        # Update values in the new 'disease' column
        adata.obs["disease"] = adata.obs["disease"].replace({"grp1": "AD", "grp2": "CTL"})
        
        # Save the modified .h5ad file (overwrite or save as new file)
        modified_file_path = os.path.join(data_dir, filename)
        adata.write_h5ad(modified_file_path, compression="gzip")
        
        print(f"Processed: {filename}")

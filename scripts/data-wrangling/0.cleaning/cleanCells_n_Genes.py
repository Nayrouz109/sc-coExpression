
#!~/anaconda3/envs/homl3/bin/python

import os
from pathlib import Path
import scanpy as sc
import sys

# Add path to the script
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/cleaning')
# Import function
from cleanCellsnGenes import clean_ctmat

# Paths
#input_path = Path("/space/scratch/nairuz-rewiring/data/0.prcsd-raw")
#output_path = Path("/space/scratch/nairuz-rewiring/data/1.prcsd-clean")
input_path = Path("/cosmos/data/project-data/NW-rewiring/data/0.prcsd-raw") #Jun 10th, 2025 
output_path = Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean")

# Ensure the output directory exists
output_path.mkdir(parents=True, exist_ok=True)

# Loop through all h5ad files in the input path
for file in input_path.glob("Mic_ROSMAP-Micro*.h5ad"):
#for file in input_path.glob("*.h5ad"):
    output_file = output_path / file.name

    # Skip processing if the file already exists in the output directory
    if output_file.exists():
        print(f"Skipping file {file.name}: already exists in the output directory.")
        continue
        
    print(f"Processing file: {file.name}")
    
    # Read the h5ad file
    adata = sc.read_h5ad(file)
    print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} genes.")
    
    # Clean the data
    cleaned_adata = clean_ctmat(adata, gene_thr=0.02, sample_thr=0.02)
    print(f"Cleaned data now has {cleaned_adata.n_obs} cells and {cleaned_adata.n_vars} genes.")
    

    
    cleaned_adata.write_h5ad(output_file, compression="gzip")
    print(f"Saved cleaned data to: {output_file}")
    print()

print("All files have been processed successfully.")


#chmod +x process_h5ad_files.py
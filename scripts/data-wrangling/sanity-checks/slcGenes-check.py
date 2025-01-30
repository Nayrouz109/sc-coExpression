
#!~/anaconda3/envs/homl3/bin/python
import sys
import anndata as ad
from pathlib import Path


input_dir = Path('/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent')
mismatched_files = []

# Initialize a variable to hold the expected number of genes (use the first file to set the expectation)
expected_num_genes = None


for file_path in input_dir.glob("*.h5ad"):
    print(f"Processing file: {file_path.name}")


    adata = ad.read_h5ad(file_path)
    
    # Check if this is the first file to set the expected number of genes
    if expected_num_genes is None:
        expected_num_genes = adata.n_vars
        print(f"Setting the expected number of genes to {expected_num_genes} based on the first file.")
    
    # Check if the current file has the expected number of genes
    if adata.n_vars != expected_num_genes:
        print(f"Warning: {file_path.name} has {adata.n_vars} genes, which does not match the expected {expected_num_genes}.")
        mismatched_files.append(file_path.name)

# Report the mismatched files
if mismatched_files:
    print("\nThe following files have a mismatched number of genes:")
    for file in mismatched_files:
        print(f"- {file}")
else:
    print("\nAll files have the same number of genes.")

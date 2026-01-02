#!~/anaconda3/envs/homl3/bin/python

import os
import sys
import json
from pathlib import Path
import anndata as ad

# Add path to the script & import functions  
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/filter-cpm')
from slcGenes_cpm import slcGenes_and_cpm_adata

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/DGE')
from create_PB import create_pseudo_bulk

# Define input and output directories
input_dir  = Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/diag-disease")
output_dir = Path('/cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/PB')
output_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

# Load gene sets
genes_sets = '/home/nelazzabi/rewiring/data/genes/genes-listsV2.json'
with open(genes_sets, "r") as fp:
    genes_set = json.load(fp)

# Step 3. Decide on the intersection set of less stringent genes
values = list(genes_set["less_stringent"].values())
intersection_set = set(values[0])
for value in values[1:]:
    intersection_set &= set(value)
intersection_set = list(intersection_set)

# Step 4. Process and save files
if __name__ == "__main__":
    for file_path in input_dir.glob("*.h5ad"):  # Loop over .h5ad files
        print(f"\nProcessing file: {file_path.name}")

        # Load AnnData object
        ann_data = ad.read_h5ad(file_path)

        # Apply filtering 
        filtered_data = slcGenes_and_cpm_adata(
            ann_data,
            intersection_set,
            cpm_normalize=False,
            log_transform=False
        )

        # Create pseudobulk data
        print("creating PB") 
        pseudobulk_expr, pseudobulk_meta = create_pseudo_bulk(filtered_data)

        # Save to CSV
        pseudobulk_expr.T.to_csv(output_dir / f"{file_path.stem}_expr.csv")
        pseudobulk_meta.to_csv(output_dir / f"{file_path.stem}_meta.csv", index=False)

    print("\nProcessing completed for all files.")












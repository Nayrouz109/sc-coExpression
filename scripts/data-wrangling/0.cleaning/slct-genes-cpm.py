#!~/anaconda3/envs/homl3/bin/python

import os
import sys
import json
from pathlib import Path
import anndata as ad

# Add path to the script & import functions  
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/filter-cpm')
from slcGenes_cpm import slcGenes_and_cpm_adata

# Define input and output directories
#input_dir = Path("/space/scratch/nairuz-rewiring/data/1.prcsd-clean")
#output_dir = Path('/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent')

#input_dir  = Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/diag-disease") #April 17th 2025 
#output_dir = Path('/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/union-lessStringent') #April 17th 2025 

#input_dir  = Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean") #June 10th, 2025 - mic study
#output_dir = Path('/cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/mic') #June 10th, 2025 - mic study

input_dir  = Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean") #Sep 22nd, 2025 - xSUbjectxType study
output_dir = Path('/cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/allCells') #Sep 22nd, 2025 - xSUbjectxType study 


output_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

# Load gene sets
genes_sets = '/home/nelazzabi/rewiring/data/genes/genes-listsV2.json'
#genes_sets = '/home/nelazzabi/rewiring/data/genes/genes-listsV3.json' #June 10th, 2025 - mic study
with open(genes_sets, "r") as fp:
    genes_set = json.load(fp)

# Step 3. Decide on the intersection set of less stringent genes
values = list(genes_set["less_stringent"].values())
intersection_set = set(values[0])
for value in values[1:]:
    intersection_set &= set(value)
intersection_set = list(intersection_set)
union_set = list(set().union(*genes_set["less_stringent"].values()))
mic_set = genes_set["less_stringent"]["Mic"]

# Step 4. Process and save files
if __name__ == "__main__":
    #for file_path in input_dir.glob("Mic*.h5ad"):  # Loop over .h5ad files
    for file_path in input_dir.glob("*.h5ad"):  # Loop over .h5ad files
        
        output_file = output_dir / file_path.name  # Retain original file name

        if output_file.exists():
            print(f"Skipping {file_path.name} — already processed.")
            continue

        print(f"\nProcessing file: {file_path.name}")

        # Load AnnData object
        #ann_data = ad.read_h5ad(file_path)

        # Apply filtering and CPM normalization
        filtered_data = slcGenes_and_cpm_adata(
            #ann_data,
            file_path,
            #union_set, 
            intersection_set,
            #mic_set, 
            cpm_normalize=False, 
            log_transform=False
        )
        
        # Save the new AnnData object
        filtered_data.write_h5ad(output_file, compression="gzip")
        print(f"Saved processed file to: {output_file}")

    print("\nProcessing completed for all files.")












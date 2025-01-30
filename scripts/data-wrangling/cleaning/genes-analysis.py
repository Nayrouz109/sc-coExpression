#!~/anaconda3/envs/homl3/bin/python

import os
import sys
import json
from pathlib import Path
import anndata as ad

# Add path to the script & import functions  
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/cleaning')
from extractGenes_perStudy import process_h5ad_files
from extractGenesIntersects import analyze_gene_dict

#input_dir = Path('/space/scratch/nairuz-rewiring/data/0.prcsd-raw')
input_dir = Path("/space/scratch/nairuz-rewiring/data/1.prcsd-clean")
output_dir = Path('/home/nelazzabi/rewiring/data/genes') 
output_file = '/home/nelazzabi/rewiring/data/genes/genes-listsV2.json'

# Utility function for JSON serialization
def convert_tuple_keys(obj):
    if isinstance(obj, dict):
        return {str(k): convert_tuple_keys(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_tuple_keys(item) for item in obj]
    else:
        return obj

# Main workflow
if __name__ == "__main__":

    # Step 1: Extract genes per study
    print("\nExtracting genes...")
    gene_dict = process_h5ad_files(input_dir)

    # Step 2: Extract intersect sets and pairs
    print("\nAnalyzing gene intersections...")
    results = analyze_gene_dict(gene_dict)
    results_serializable = convert_tuple_keys(results)

    # Step 3: Save results
    print("\nSaving results to JSON...")
    with open(output_file, 'w') as f:
        json.dump(results_serializable, f, indent=4)

    # Step 4: Print summary of results
    print("\nSummary of stringent genes per cell type:")
    for cell_type, genes in results["stringent"].items():
        print(f"{cell_type}: {len(genes)} genes")

    print("\nSummary of less stringent genes per cell type:")
    for cell_type, genes in results["less_stringent"].items():
        print(f"{cell_type}: {len(genes)} genes")

    print("\nPairwise intersections of less stringent genes:")
    for pair, intersect_genes in results["pairwise_intersections"].items():
        print(f"{pair}: {len(intersect_genes)} genes")

    print("\nProcessing completed successfully.")

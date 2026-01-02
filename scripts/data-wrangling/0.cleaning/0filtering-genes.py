#!~/anaconda3/envs/homl3/bin/python

import os
import numpy as np
import pandas as pd
import anndata as ad
import json


input_dir = '/space/scratch/nairuz-rewiring/data/0.prcsd-raw'
output_dir = '/home/nelazzabi/rewiring/data/genes'

# Make sure output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

gene_dict = {}


def extract_metadata_from_filename(filename):
    # Assuming the filename format: "CellType_StudyName_Year.h5ad"
    name_parts = filename.split('_')
    cell_type = name_parts[0] 
    study_name = name_parts[-1].split('.')[0]  
    return cell_type, study_name



def process_h5ad_files(directory):

    for filename in os.listdir(directory):
        if filename.endswith('.h5ad'):
            file_path = os.path.join(directory, filename)
            print(f"\nProcessing file: {filename}")

            try:
                adata = ad.read_h5ad(file_path)
            except Exception as e:
                print(f"Error reading {filename}: {e}")
                continue

            # Print the number of genes before filtering
            total_genes = adata.var_names.size
            print(f"Total genes in {filename}: {total_genes}")

            # Filter out genes with zero expression in all cells
            non_zero_genes = np.asarray(adata.X.sum(axis=0)) > 0
            non_zero_genes = non_zero_genes.flatten()  # Flatten to 1D array
            valid_genes = adata.var_names[non_zero_genes]

            # Print the number of genes after filtering
            filtered_genes_count = valid_genes.size
            print(f"Valid genes after filtering in {filename}: {filtered_genes_count}")

            # Extract metadata (cell type, study)
            cell_type, study_name = extract_metadata_from_filename(filename)

            # Add the valid genes to the dictionary under cell type -> study name
            if cell_type not in gene_dict:
                gene_dict[cell_type] = {}  # Add cell type if not already present
            if study_name not in gene_dict[cell_type]:
                gene_dict[cell_type][study_name] = []  # Add study name if not already present

            # Add valid genes to the dictionary
            if filtered_genes_count > 0:
                gene_dict[cell_type][study_name] = valid_genes.tolist()
            else:
                print(f"No valid genes found in {filename}")

    return gene_dict

# Main execution
if __name__ == '__main__':
    print("Starting gene filtering process...\n")

    # Process .h5ad files in the given directory
    gene_dict = process_h5ad_files(input_dir)

    # Save the resulting genes to a CSV file
    output_file = os.path.join(output_dir, 'filtered_genes.json')

    print("\nSaving results...")
    with open(output_file, 'w') as f:
        json.dump(gene_dict, f)
        print(f"Gene dictionary saved to {output_file}")


    print(f"Filtering complete. Results saved to {output_file}\n")
    print("Process complete.")










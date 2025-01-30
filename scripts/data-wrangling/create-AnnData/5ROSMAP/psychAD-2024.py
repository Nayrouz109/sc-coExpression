

#!~/anaconda3/envs/homl3/bin/python
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import csr_matrix


# Define input and output directories =============================
input_dir  = '/cosmos/data/downloaded-data/AD-meta/psychAD-2024/data/'
output_dir = '/cosmos/data/downloaded-data/AD-meta/psychAD-2024/data/annData-cellType/'

print("loading AnnData from:", input_dir)
data = ad.read_h5ad(os.path.join(input_dir, "psychAD_snRNAseq_rawCounts.h5ad"))
data.obs = data.obs.rename(columns={"Individual_ID": "patientID"})


cell_types = {
    "Exc": data.obs['class'].isin(['EN']),
    "Inh": data.obs['class'].isin(['IN']),
    "Ast": data.obs['class'].isin(['Astro']),
    "Oli": data.obs['class'].isin(['Oligo']),
    "Opc": data.obs['class'].isin(['OPC']),
    "Mic": data.obs['subclass'].isin(['Micro']),
}
sources = data.obs["Source"].unique().tolist()

# Iterate over cell types and create subsets =======================
for cell_type, condition in cell_types.items():
    for source in sources:
        combined_condition = condition & (data.obs["Source"] == source)
        print(f"Processing cell type: {cell_type}, Source: {source}")

        subset_data = data[combined_condition]
        reduced_data = ad.AnnData(X=csr_matrix(subset_data.X), obs=subset_data.obs, var=data.var)
        reduced_data.obs = reduced_data.obs.astype(str) 

        # Format the file name for saving
        file_name = f"{cell_type}_{source}-2024.h5ad"  
        file_path = output_dir + file_name  # Create the full file path
    
        # Save the reduced AnnData object: Apply Compression When Saving:
        reduced_data.write_h5ad(file_path, compression="gzip")
        print(f"Saved subset for {cell_type} to {file_path}")

# Optionally, you can delete variables to free up memory
del data, subset_data, reduced_data


print("Processing complete.")















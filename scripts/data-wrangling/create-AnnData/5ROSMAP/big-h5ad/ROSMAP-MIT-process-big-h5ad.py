



#!~/anaconda3/envs/homl3/bin/python
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import csr_matrix


# Define input and output directories =============================
input_dir  = '/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT/data/0-source/expr/h5ad/'
output_dir = '/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT/data/0-source/expr/h5ad/annData-cellType/'

print("loading AnnData from:", input_dir)
data = ad.read_h5ad(os.path.join(input_dir, "PFC427_raw_data.h5ad"))
data.obs = data.obs.rename(columns={"individual_ID": "patientID"})



# Define the cell types and conditions for subsetting ===================
cell_types = {
    "Exc": data.obs['major_cell_type'].isin(['Exc']),
    "Inh": data.obs['major_cell_type'].isin(['Inh']),
    "Ast": data.obs['major_cell_type'].isin(['Ast']),
    "Oli": data.obs['major_cell_type'].isin(['Oli']),
    "Opc": data.obs['major_cell_type'].isin(['OPC']),
    "Mic": data.obs['major_cell_type'].isin(['Mic']),
}


# Iterate over cell types and create subsets =======================
for cell_type, condition in cell_types.items():
    print(f"Processing cell type: {cell_type}")

    # Subset the AnnData object without storing extra variables
    subset_data = data[condition]  # No copy to save memory
    reduced_data = ad.AnnData(X=csr_matrix(subset_data.X), obs=subset_data.obs, var=data.var)
    reduced_data.obs = reduced_data.obs.astype(str)

    # Format the file name for saving
    file_name = f"{cell_type}_ROSMAP-Mathys-2023.h5ad"  
    file_path = output_dir + file_name  # Create the full file path
    
    # Save the reduced AnnData object: Apply Compression When Saving:
    reduced_data.write_h5ad(file_path, compression="gzip")
    print(f"Saved subset for {cell_type} to {file_path}")

# Optionally, you can delete variables to free up memory
del data, subset_data, reduced_data


print("Processing complete.")

#chmod +x process_ann_data.py
# python process_ann_data.py




#!~/anaconda3/envs/homl3/bin/python
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# Print the version of the anndata package
print("AnnData version:", ad.__version__)

# Define input and output directories
input_dir = '/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/expr/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad'
output_dir = '/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/expr/annData-cellTpe/'

# Load the AnnData object
print("Loading AnnData from:", input_dir)
data = ad.read_h5ad(input_dir)
print("Loaded AnnData with shape:", data.shape)

# Define the cell types and conditions for subsetting
cell_types = {
    "Exc": data.obs['Class'].isin(['Neuronal: Glutamatergic']),
    "Inh": data.obs['Class'].isin(['Neuronal: GABAergic']),
    "Ast": data.obs['Subclass'].isin(['Astrocyte']),
    "Oli": data.obs['Subclass'].isin(['Oligodendrocyte']),
    "Opc": data.obs['Subclass'].isin(['OPC']),
    "Mic": data.obs['Subclass'].isin(['Microglia-PVM']),
}

# Iterate over cell types and create subsets
for cell_type, condition in cell_types.items():
    print(f"Processing cell type: {cell_type}")

    # Subset the AnnData object without storing extra variables
    subset_data = data[condition]  # No copy to save memory
    reduced_data = ad.AnnData(X=csr_matrix(subset_data.X), obs=subset_data.obs, var=data.var)

    # Subset the AnnData object without storing extra variables
    #reduced_data.layers.clear()
    #reduced_data.obsm.clear()
    #reduced_data.varm.clear()
    #reduced_data.uns.clear()

    # Format the file name for saving
    file_name = f"{cell_type}_SEA-AD-MTG-2024.h5ad"  
    file_path = output_dir + file_name  # Create the full file path
    
    # Save the reduced AnnData object: Apply Compression When Saving:
    #reduced_data.write_h5ad(file_path)  
    reduced_data.write_h5ad(file_path, compression="gzip")
    print(f"Saved subset for {cell_type} to {file_path}")

# Optionally, you can delete variables to free up memory
del data, subset_data, reduced_data

print("Processing complete.")

#chmod +x process_ann_data.py
# python process_ann_data.py







#!~/anaconda3/envs/homl3/bin/python
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import csr_matrix


# Define input and output directories =============================
input_dir = '/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT-multiregion-2024/data/'
output_dir = '/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT-multiregion-2024/data/annData-cellType/'

print("Loading AnnData from:", input_dir)
data = ad.read_h5ad(os.path.join(input_dir, "all_brain_regions_filt_preprocessed_scanpy_fullmatrix.h5ad"))
cells_meta = pd.read_csv(os.path.join(input_dir, "all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv"), delimiter='\t')


# sanity check ===================================================
missing_barcodes = set(cells_meta['barcode'].values) - set(data.obs_names)
len(missing_barcodes) #should be 0 

# filter to 1.3 Mill cells ===================================================
desired_cells = list(set(cells_meta['barcode'].values))
data_subset = data[data.obs_names.isin(desired_cells)] 
del data 

# merge cell-level meta data ===================================================
print("merging cell-level meta data")
# 1. Set 'barcode' as the index in cells_meta to prepare for merging
cells_meta.set_index('barcode', inplace=True)
# 2. Make sure that data_subset.obs is indexed by the barcodes (if it isn't already)
data_subset.obs.set_index(data_subset.obs_names, inplace=True)
# 3. Merge data_subset.obs with cells_meta, using the index (barcode) for alignment
data_subset.obs = data_subset.obs.join(cells_meta, how="left")

data_subset.obs = data_subset.obs.rename(columns={"projid": "patientID"})



# Define the cell types and conditions for subsetting ===================
cell_types = {
    "Exc": data_subset.obs['minor.celltype'].isin(['Exc']),
    "Inh": data_subset.obs['minor.celltype'].isin(['Inh']),
    "Ast": data_subset.obs['minor.celltype'].isin(['Ast']),
    "Oli": data_subset.obs['minor.celltype'].isin(['Oli']),
    "Opc": data_subset.obs['minor.celltype'].isin(['Opc']),
    "Mic": data_subset.obs['minor.celltype'].isin(['Mic']),
}


# Iterate over cell types and create subsets =======================
for cell_type, condition in cell_types.items():
    print(f"Processing cell type: {cell_type}")

    # Subset the AnnData object without storing extra variables
    subset_data = data_subset[condition]  # No copy to save memory
    reduced_data = ad.AnnData(X=csr_matrix(subset_data.X), obs=subset_data.obs, var=data_subset.var)
    reduced_data.obs = reduced_data.obs.astype(str)

    # Format the file name for saving
    file_name = f"{cell_type}_ROSMAP-Mathys-2024.h5ad"  
    file_path = output_dir + file_name  # Create the full file path
    
    # Save the reduced AnnData object: Apply Compression When Saving:
    reduced_data.write_h5ad(file_path, compression="gzip")
    print(f"Saved subset for {cell_type} to {file_path}")

# Optionally, you can delete variables to free up memory
del data_subset, subset_data, reduced_data


print("Processing complete.")

#chmod +x process_ann_data.py
# python process_ann_data.py
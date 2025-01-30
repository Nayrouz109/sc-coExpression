

#!~/anaconda3/envs/homl3/bin/python
import numpy as np
import pandas as pd
import anndata as ad
import os
import scipy.sparse as sp
import scipy.io
from scipy.sparse import csr_matrix


# read in data =============================
input_dir  = '/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/'
output_dir = '/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/annData-cellType/'

print("Loading data")
counts = scipy.io.mmread(os.path.join(input_dir, "ROSMAP_Brain.snRNAseq_counts_sparse_format_20230420.csv")).tocsc()
genes  = pd.read_csv(os.path.join(input_dir, "ROSMAP_Brain.snRNAseq_metadata_genes_20230420.csv"))
cells_meta = pd.read_csv(os.path.join(input_dir, "ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv"))

# build annData =============================
print("building AnnData")
counts = counts.transpose()  # make it cells x genes
genes = genes[["x"]]
genes.columns = ["gene_names"]

adata = ad.AnnData(X=counts, obs=cells_meta, var=genes)
adata.obs_names = cells_meta["cell_name"]
adata.var_names = genes["gene_names"]

adata.obs = adata.obs.rename(columns={"specimenID": "patientID"})

# Define the cell types and conditions for subsetting ===================
cell_types = {
    "Exc": adata.obs['Cell Type'].isin(['Exc']),
    "Inh": adata.obs['Cell Type'].isin(['Inh']),
    "Ast": adata.obs['Cell Type'].isin(['astrocytes']),
    "Oli": adata.obs['Cell Type'].isin(['oligodendrocytes']),
    "Opc": adata.obs['Cell Type'].isin(['opcs']),
    "Mic": adata.obs['Cell Type'].isin(['microglia']),
}

# Iterate over cell types and create subsets =======================
for cell_type, condition in cell_types.items():
    print(f"Processing cell type: {cell_type}")

    # Subset the AnnData object without storing extra variables
    subset_data = adata[condition]  # No copy to save memory
    reduced_data = ad.AnnData(X=csr_matrix(subset_data.X), obs=subset_data.obs, var=adata.var)
    reduced_data.obs = reduced_data.obs.astype(str)

    # Format the file name for saving
    file_name = f"{cell_type}_ROSMAP-DeJager-2023.h5ad"  
    file_path = output_dir + file_name  # Create the full file path
    
    # Save the reduced AnnData object: Apply Compression When Saving:
    reduced_data.write_h5ad(file_path, compression="gzip")
    print(f"Saved subset for {cell_type} to {file_path}")

# Optionally, you can delete variables to free up memory
del adata, subset_data, reduced_data


print("Processing complete.")




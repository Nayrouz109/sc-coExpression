

#!~/anaconda3/envs/homl3/bin/python
import pandas as pd
import scipy.io
import anndata as ad
import os
import scipy.io
from scipy.sparse import csr_matrix

# Paths to the saved data
path = "/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components"


### create annData ============================================================
# Load data 
print(f"Reading count matrix") 
counts = scipy.io.mmread(os.path.join(path, "countsV2.mtx")).tocsc()
genes = pd.read_csv(os.path.join(path, "genes.csv"))
cells = pd.read_csv(os.path.join(path, "cell_ids.csv"))
obs = pd.read_csv(os.path.join(path, "meta_data.csv"), index_col=0)

# Ensure row and column names match
counts = counts.transpose()  # make it cells x genes
genes.columns = ["gene_names"]
cells.columns = ["cell_ids"]

# Create AnnData object
adata = ad.AnnData(X=counts, obs=obs, var=genes)
adata.obs_names = cells["cell_ids"]
adata.var_names = genes["gene_names"]


### create annData cell subsets ============================================================
output_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw/"
cell_types = {
    "Exc": adata.obs['cell_type'].isin(['excitatory']),
    "Inh": adata.obs['cell_type'].isin(['inhibitory']),
    "Ast": adata.obs['cell_type'].isin(['astrocyte']),
    "Oli": adata.obs['cell_type'].isin(['oligodendrocyte']),
    "Opc": adata.obs['cell_type'].isin(['opc']),
    "Mic": adata.obs['cell_type'].isin(['microglia']),
}

# Iterate over cell types and create subsets
for cell_type, condition in cell_types.items():
    print(f"Processing cell type: {cell_type}")

    # Subset the AnnData object without storing extra variables
    subset_data = adata[condition]  # No copy to save memory
    reduced_data = ad.AnnData(X=csr_matrix(subset_data.X), obs=subset_data.obs, var=adata.var)

    # Sanity Check 1: Dimension Consistency
    assert reduced_data.X.shape[0] == reduced_data.obs.shape[0], "Mismatch between cells in X and obs."
    assert reduced_data.X.shape[1] == reduced_data.var.shape[0], "Mismatch between genes in X and var."
    print("Dimension check passed.")

    # Sanity Check 2: Name Alignment
    assert (reduced_data.obs_names == reduced_data.obs.index).all(), "Mismatch between obs_names and obs index."
    assert (reduced_data.var_names == reduced_data.var.index).all(), "Mismatch between var_names and var index."
    print("Name alignment check passed.")

    # Sanity Check 3: Matrix Sparsity
    if csr_matrix(subset_data.X).nnz > 0:
        print(f"{cell_type} data is sparse, with shape:", reduced_data.X.shape)
    else:
        print(f"{cell_type} data is dense, consider using a sparse format.")

    # Sanity Check 5: Expression Matrix Range Check
    min_val = reduced_data.X.min()
    max_val = reduced_data.X.max()
    print(f"Expression matrix values range from {min_val} to {max_val}")
    assert min_val >= 0, "Expression matrix contains negative values."

    # Format the file name for saving
    file_name = f"{cell_type}_Rexach-2024.h5ad"  
    file_path = output_dir + file_name  # Create the full file path
    
    # Save the reduced AnnData object: Apply Compression When Saving:
    reduced_data.write_h5ad(file_path, compression="gzip")
    print(f"Saved subset for {cell_type} to {file_path}")

# Optionally, you can delete variables to free up memory
del adata, subset_data, reduced_data

print("Processing complete.")


#chmod +x <script>.py
# python <script>.py

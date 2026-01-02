

# %% [markdown]
# ## data import and setup 
import numpy as np
import pandas as pd 
import anndata as ad 
from scipy.sparse import csr_matrix
print(ad.__version__)



data = '/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad'
data = ad.read_h5ad(data)

data.X     #data matrix, numpy.ndarray or scipy.sparse 
data.obs   # pandas.dataframe 


#%%
cell_types = {
    "Exc": data.obs['Class'].isin(['Neuronal: Glutamatergic']),
    "Inh": data.obs['Class'].isin(['Neuronal: GABAergic']),
    "Ast": data.obs['Subclass'].isin(['Astrocyte']),
    "Oli": data.obs['Subclass'].isin(['Oligodendrocyte']),
    "Opc": data.obs['Subclass'].isin(['OPC']),
    "Mic": data.obs['Subclass'].isin(['Microglia-PVM']),
}

output_dir = '/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/annData-cellTpe/'

cell_type = next(iter(cell_types))

reduced_data = anndata.AnnData(X=subset_data.X, obs=subset_data.obs, var=subset_data.var)

# %%
for cell_type, condition in cell_types.items():
    subset_data = data[condition].copy()  # Create a copy of the subsetted AnnData object
    file_name = f"{cell_type}_SEA-AD-DLPFC-2024.h5ad"  # Format the file name
    file_path = output_dir + file_name  # Create the full file path
    subset_data.write_h5ad(file_path)  # Save the subsetted AnnData object
    print(f"Saved subset for {cell_type} to {file_path}"
)
    


# %%    
mic = "/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/expr/annData-cellTpe/Mic_SEA-AD-DLPFC-2024.h5ad"
data = ad.read_h5ad(mic)
# %%
micDLPFC = "/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/annData-cellTpe/Mic_SEA-AD-DLPFC-2024.h5ad"
data = ad.read_h5ad(micDLPFC)

data.obs.columns
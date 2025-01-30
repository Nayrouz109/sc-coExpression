





#### testing inner fuctions ==========================
### packages 
import sys
import os
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
from joblib import Parallel, delayed
import bottleneck 
from pathlib import Path
import multiprocessing
import scipy.sparse
 




 





# Add the path to source the function
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV0_disease import compute_gene_correlation
from corOnly_parallelV4_disease import compute_gene_correlation
from corOnly_parallelV0_disease import sparse_corr
from corOnly_parallelV0_disease import rank
from saveCorDict import save_to_hdf5


### arguments 
input_dir = "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease"
output_dir = "/space/scratch/nairuz-rewiring/coExpr/AD-corMtx"

correlation_typeT =  "pearson" 
replace_nansT = True
analysis_typeT =  "disease"  
categoriesT = ['AD', 'CTL']
file_pathT = "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/Mic_ROSMAP-Mathys-2023.h5ad"



adata = sc.read_h5ad("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/Mic_ROSMAP-Mathys-2023.h5ad")
meta = adata.obs

# Randomly sample 10 unique patient IDs for each disease category
sampled_patients = meta.groupby('disease')['patientID'].apply(
    lambda x: np.random.choice(x.unique(), 5, replace=False)
).explode()

# Subset the AnnData object
adata_subset = adata[meta['patientID'].isin(sampled_patients)].copy()




cor = compute_gene_correlation(
                adata=adata_subset,
                correlation_type= correlation_typeT,
                replace_nans= replace_nansT,
                min_cells_threshold=20,
                max_workers=6,
                analysis_type= analysis_typeT,
                categories= categoriesT,
                file_path=file_pathT
)


#### testing outer fuctions ==========================
python run_compute_cor-disease.py \
    --input_dir /space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease \
    --output_dir /space/scratch/nairuz-rewiring/coExpr/AD-corMtx \
    --correlation_type pearson \
    --replace_nans True \
    --analysis_type disease \
    --categories AD CTL


python run_compute_cor-disease.py \
    --input_dir /space/scratch/nairuz-rewiring/data/testing \
    --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease \
    --correlation_type pearson \
    --replace_nans True \
    --analysis_type disease \
    --categories AD CTL



### notes to edit the compute_gene_correlation()
# 1. make sure the script works for both disease and non-disease
# 2. to make things easy to follow, create a new function that aggregates and summarize results 


#### testing rank function ==========================
from scipy.sparse import csr_matrix
subset = adata[adata.obs['disease'] == "AD"]
patient_matrix = subset[subset.obs['patientID'] == "11336574"].X


adata_subset = adata[adata.obs["patientID"] == "41285665"]
adata_subset.shape
subset_data = adata_subset.X
subset_data = csr_matrix(subset_data)



corr_matrix = sparse_corr(subset_data)
np.fill_diagonal(corr_matrix, 1)
rank_normalized_matrix = rank(corr_matrix)
rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = bottleneck.nanmean(rank_normalized_matrix)




### check the output 
edgeAD = pd.read_csv("/cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list/Mic_Rexach-2024_AD_corEdgeLists.csv", index_col=[0, 1])
edgeCT = pd.read_csv('/cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list/Mic_Rexach-2024_CTL_corEdgeLists.csv', index_col=[0, 1])
aggAD = ad.read_h5ad("/cosmos/data/project-data/NW-rewiring/coExpr/disease/agg-mtx/Mic_Rexach-2024_AD_corAggSparse.h5ad")
aggCT = ad.read_h5ad("/cosmos/data/project-data/NW-rewiring/coExpr/disease/agg-mtx/Mic_Rexach-2024_CTL_corAggSparse.h5ad")


input_dir = Path("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease")


for file_path in input_dir.glob("Inh*.h5ad"):
    print(file_path) 
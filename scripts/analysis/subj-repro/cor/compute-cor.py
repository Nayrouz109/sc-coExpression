#!~/anaconda3/envs/homl3/bin/python
import os
import pandas as pd
import numpy as np
import anndata as ad 
import sys
import scipy.sparse
from scipy.stats import spearmanr
import random
from anndata import read_h5ad
from joblib import Parallel, delayed

sys.path.append('/home/nelazzabi/rewiring/scripts/analysis/subj-repro')
from functions import process_and_save


# Paths for data
input_files = {
    "2023": "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag/Exc_ROSMAP-Mathys-2023.h5ad",
    "2024": "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag/Exc_ROSMAP-Mathys-2024.h5ad", 
    "random": "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag/Opc_ling-2024.h5ad"
}
patient_ids_path = "/home/nelazzabi/rewiring/scripts/analysis/subj-repro/patientIDs.csv"
output_dirs = {
    "raw": "/space/scratch/nairuz-rewiring/coExpr/subj-repro/raw-cor",
    "rank": "/space/scratch/nairuz-rewiring/coExpr/subj-repro/rank-norm-raw"
}


# Read patient IDs
df_patient_ids = pd.read_csv(patient_ids_path)
patient_ids = df_patient_ids['patientID'].astype(str).tolist()

# Process dataset
#results_2023 = process_and_save(input_files["2023"], patient_ids, "ROSMAP-Mathys-2023",output_dirs )
#results_2024 = process_and_save(input_files["2024"], patient_ids, "ROSMAP-Mathys-2024", output_dirs, brain_region="PFC" ) 
results_rndm  = process_and_save(input_files["random"], 
                                 random.sample(read_h5ad(input_files["random"]).obs["patientID"].unique().astype(str).tolist(), 20), 
                                 "ling-2024",
                                 output_dirs )

print("Processing complete. Results saved to specified directories.\n")

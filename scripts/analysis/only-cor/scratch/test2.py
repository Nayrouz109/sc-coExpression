


#!~/anaconda3/envs/homl3/bin/python
import sys
import os
from pathlib import Path
import argparse
import scanpy as sc
import numpy as np
import gc 
import psutil

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV9_disease import compute_gene_correlation

input_dir = Path("/space/scratch/nairuz-rewiring/data/testing/test")
output_dir = "/space/scratch/nairuz-rewiring/coExpr/AD-corMtx"
categoriesT = ['AD', 'CTL']
correlation_typeT =  "pearson" 
replace_nansT = True
analysis_typeT =  "disease"  

for file_path in input_dir.glob("*.h5ad"):
    print(f"\nProcessing file: {file_path}")

    try:
        adata = sc.read_h5ad(file_path)
        meta = adata.obs
        sampled_patients = meta.groupby('disease')['patientID'].apply(lambda x: np.random.choice(x.unique(), 7, replace=False)).explode()

        # Subset the AnnData object
        adata = adata[meta['patientID'].isin(sampled_patients)].copy()
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        continue

    try:
        for category in categoriesT: 
            compute_gene_correlation(
                adata=adata,
                correlation_type=correlation_typeT,
                replace_nans=replace_nansT,
                min_cells_threshold=20,
                category=category,
                file_path=file_path,
                output_dir=output_dir
            )

    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        continue
    del adata
    gc.collect()





























import pandas as pd
import scipy.io
import anndata as ad
import os

adata = ad.read_h5ad('Mic_ROSMAP-Micro.h5ad')
obs = adata.obs

obs = obs.rename(columns = {'patientID' : 'patientID_original'})

obs['patientID'] = obs['patientID_original'].astype(str) + '_' + obs['brainRegion'].astype(str)

obs['disease'] = obs['ADdiag3types'].map({
    'earlyAD': 'AD',
    'lateAD': 'AD',
    'nonAD': 'CTL'
}).fillna('NA')

obs['diagnosis'] = obs['disease'].map({
    'AD': 'yes',
    'CTL': 'no'
}).fillna('NA')


adata.obs = obs

output_path = os.path.join('/space/scratch/bxu_project/ROSMAP-Micro/Mic_ROSMAP-Micro-V1.h5ad')
adata.write_h5ad(output_path, compression="gzip")
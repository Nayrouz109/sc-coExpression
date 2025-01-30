#!~/anaconda3/envs/homl3/bin/python

import pandas as pd 
import anndata as ad 
import numpy as np
import sys

sys.path.append('/home/nelazzabi/rewiring/scripts/analysis/subj-repro')
from fnct_compute_OR import perform_comparisons

#### create the grid  =========================================
summary_2024 = pd.read_csv("/space/scratch/nairuz-rewiring/coExpr/subj-repro/raw-cor/summary_ROSMAP-Mathys-2024_PFC.csv")
summary_2023 = pd.read_csv("/space/scratch/nairuz-rewiring/coExpr/subj-repro/raw-cor/summary_ROSMAP-Mathys-2023_NA.csv")
grid = pd.concat([summary_2024, 
                  summary_2023], axis=0) #axis = 0 row-wise 
grid.reset_index(drop=True, inplace=True)
grid = grid.sort_values(by='PatientID', ascending=True, na_position='last')

### Perform comparisons ============================================================================
results = perform_comparisons(grid)

results.to_csv('/home/nelazzabi/rewiring/scripts/analysis/subj-repro/results.csv', index=False)

#!~/anaconda3/envs/homl3/bin/python

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc 
import sys
from pathlib import Path


sys.path.append('/home/nelazzabi/rewiring/scripts/functions/DGE')
from create_PB import create_pseudo_bulk


# Process all files
#data_dir = Path("/space/scratch/nairuz-rewiring/data/1.prcsd-clean")
#output_dir = Path("/space/scratch/nairuz-rewiring/data/1.prcsd-clean/pseudobulk")
#files = list(data_dir.glob("*ling-2024.h5ad"))

#data_dir =   Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/diag-disease") #feb 18th 
#files = list(data_dir.glob("*.h5ad"))
#output_dir = Path("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk")

data_dir =   Path("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease") #march 18th 
files = list(data_dir.glob("*.h5ad"))
output_dir = Path("/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/pseudobulk/raw")

output_dir.mkdir(exist_ok=True)

for file in files:
    adata = ad.read_h5ad(file)

    # Create pseudobulk data
    print(f"processing {file.stem}")
    pseudobulk_expr, pseudobulk_meta = create_pseudo_bulk(adata)

    # Save to CSV
    pseudobulk_expr.T.to_csv(output_dir / f"{file.stem}_expr.csv")
    pseudobulk_meta.to_csv(output_dir / f"{file.stem}_meta.csv", index=False)

    print(f"Processed and saved pseudobulk data for {file.stem}\n")




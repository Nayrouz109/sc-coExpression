#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

"""
Executable script to generate link tables and compute Fisher's exact tests
for reproducibility across cell types.
"""

import os
import argparse
import pandas as pd

import sys
sys.path.append("/home/nelazzabi/rewiring/scripts/functions/xCellxSubject")
from lnk_fisherV2 import process_lnktbles, compute_fisher_XcellXsubj

def main():
    parser = argparse.ArgumentParser(description="Generate link tables and compute Fisher's test")
    parser.add_argument("--input_dir", required=True, help="Directory containing h5ad files")
    parser.add_argument("--lnk_output_dir", required=True, help="Directory to save link tables")
    parser.add_argument("--fisher_output_dir", required=True, help="Directory to save Fisher test results")
    parser.add_argument("--cell_types", nargs="+", default=["Mic", "Oli", "Opc", "Exc", "Inh", "Ast"], help="List of cell types to process")
    parser.add_argument("--condition", type=str, default=None, help="Condition string to filter files (e.g., 'CTL', 'AD'). If not provided, all conditions are used.")
    #parser.add_argument("--n_jobs", type=int, default=4, help="Number of parallel jobs for Fisher test")
    
    args = parser.parse_args()
    if args.cell_types == ["None"]:
        args.cell_types = None 


    # Step 1: Generate link tables
    print("=== Step 1: Generating link tables ===")
    process_lnktbles(args.input_dir, args.lnk_output_dir, args.cell_types, args.condition)
    
    # Step 2: Compute Fisher's exact tests
    print("\n=== Step 2: Computing Fisher's exact tests ===")
    compute_fisher_XcellXsubj(args.lnk_output_dir, 
                              args.fisher_output_dir, 
                              args.cell_types, 
                              args.condition)


    print("\n✅ Workflow complete!")

if __name__ == "__main__":
    main()



# Sep 30th, 2025 
##############################################################################################################
############################################ xSubjectxType ################################################### 
#python run_xCellxSubject.py \
#    --input_dir '/cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType/agg-mtx' \
#    --lnk_output_dir    '/cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType/rpr/lnkTables' \
#    --fisher_output_dir "/cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType/rpr/fisher" \
#    --cell_types None \
#    --condition CTL

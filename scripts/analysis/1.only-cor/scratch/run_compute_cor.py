#!~/anaconda3/envs/homl3/bin/python

import sys
import os
from pathlib import Path
import argparse
import scanpy as sc


# Add the path to source the function
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV0 import compute_gene_correlation
from saveCorDict import save_to_hdf5

# input and output directories 
#input_dir = Path("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag")
#output_dir = Path("/space/scratch/nairuz-rewiring/coExpr")
#output_dir = Path("/space/scratch/nairuz-rewiring/coExpr/AD-corMtx")

# prsd arguments 
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run compute_gene_correlation on a single .h5ad file.")
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Path to the input directory containing .h5ad files.")
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Path to the output directory for saving results.")
    parser.add_argument("--correlation_type", type=str, default="pearson", choices=["pearson", "spearman"],
                        help="Type of correlation to compute (default: 'pearson').")
    parser.add_argument("--replace_nans", type=bool, default=True,
                        help="Whether to replace NaNs in the correlation matrix with the median value (default: True).")
    parser.add_argument("--analysis_type", type=str, default=None, choices=["baseline", None],
                        help="Type of analysis. If 'baseline', subset data where diagnosis == 'no'.")
    args = parser.parse_args()

    for file_path in input_dir.glob("*.h5ad"):
        print(f"Processing file: {file_path}")

        try:
            adata = sc.read_h5ad(file_path)
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            continue

        try:
            cor = compute_gene_correlation(
                adata=adata,
                correlation_type=args.correlation_type,
                replace_nans=args.replace_nans,
                min_cells_threshold=20,
                max_workers=6,
                analysis_type=args.analysis_type,
                file_path=file_path
            )

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            continue

        # Save the result
        try:
            save_to_hdf5(cor, output_dir, file_path)
        except Exception as e:
            print(f"Error saving results for file {file_path}: {e}")

if __name__ == "__main__":
    main()



### example of how to run the code 
#~/anaconda3/envs/homl3/bin/python your_script_name.py \
    #--input_dir /space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag \
    #--output_dir /space/scratch/nairuz-rewiring/coExpr/AD-corMtx \
    #--correlation_type pearson \
    #--replace_nans True \
    #--analysis_type None









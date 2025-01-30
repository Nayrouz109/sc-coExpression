#!~/anaconda3/envs/homl3/bin/python
import sys
import os
from pathlib import Path
import argparse
import scanpy as sc
import gc 

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV6_disease import compute_gene_correlation


def main():
    parser = argparse.ArgumentParser(description="Run compute_gene_correlation on a single .h5ad file.")
    parser.add_argument("--input_dir", type=str, required=True, 
                        help="Path to the input directory containing .h5ad files.")
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Path to the output directory for saving results.")
    parser.add_argument("--correlation_type", type=str, default="pearson", choices=["pearson", "spearman"],
                        help="Type of correlation to compute (default: 'pearson').")
    parser.add_argument("--replace_nans", action='store_true', 
                        help="Whether to replace NaNs in the correlation matrix with the median value (default: True).")
    parser.add_argument("--analysis_type", type=str, default=None, choices=["baseline", "disease"],
                        help="Type of analysis. 'baseline' for control data, 'disease' for specified categories.")
    parser.add_argument("--categories", type=str, nargs='+', default=['AD', 'CTL'],
                        help="Categories to subset the data when analysis_type is 'disease' (default: ['AD', 'CTL']).")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    for file_path in input_dir.glob("*.h5ad"):
        print(f"\nProcessing file: {file_path}")

        try:
            adata = sc.read_h5ad(file_path)
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            continue

        if args.analysis_type == "disease":
            if not set(args.categories).issubset(adata.obs["disease"].unique()):
                print(f"File {file_path} does not contain all required categories {args.categories}. Skipping.")
                continue

        try:
            for category in args.categories: 
                compute_gene_correlation(
                    adata=adata,
                    correlation_type=args.correlation_type,
                    replace_nans=args.replace_nans,
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


if __name__ == "__main__":
    main()


# Jan 29th 
# python run_compute_cor-diseaseV0.py --input_dir /space/scratch/nairuz-rewiring/data/testing/test --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease --correlation_type pearson --replace_nans --analysis_type disease --categories AD CTL


#python run_compute_cor-diseaseV0.py \
    #--input_dir /space/scratch/nairuz-rewiring/data/testing/test \
    #--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease \
    #--correlation_type pearson \ 
    #--replace_nans \
    #--analysis_type disease \ 
    #--categories AD CTL


#if replace_nas is false 
#python run_correlation.py 
# --input_dir /path/to/input_directory 
# --output_dir /path/to/output_directory 
# --correlation_type spearman 
# --analysis_type disease 
# --categories AD CTL

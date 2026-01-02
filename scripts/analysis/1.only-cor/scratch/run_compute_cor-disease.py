#!~/anaconda3/envs/homl3/bin/python

### readMe =======================================================
### this is to compute co-expression networks in disease context 

import sys
import os
from pathlib import Path
import argparse
import scanpy as sc

# Add the path to source the function
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV5_disease import compute_gene_correlation
from saveCorDict import save_to_hdf5

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
    parser.add_argument("--analysis_type", type=str, default=None, choices=["baseline", "disease"],
                        help="Type of analysis. 'baseline' for control data, 'disease' for specified categories.")
    parser.add_argument("--categories", type=str, nargs='+', default=['AD', 'CTL'],
                        help="Categories to subset the data when analysis_type is 'disease' (default: ['AD', 'CTL']).")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

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
            cor = compute_gene_correlation(
                adata=adata,
                correlation_type=args.correlation_type,
                replace_nans=args.replace_nans,
                min_cells_threshold=20,
                #max_workers=6,
                analysis_type=args.analysis_type,
                categories=args.categories,
                file_path=file_path
            )

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")
            continue

        # Save the results for each category in the dictionary
        for category in cor:
            category_data = cor[category]
            category_file_path = Path(category_data["file_path"])
            try:
                print("saving dictionary") 
                save_to_hdf5(category_data, output_dir, category_file_path)
                print(f"Results for {category} saved successfully.")
            except Exception as e:
                print(f"Error saving results for {category} in file {file_path}: {e}")

if __name__ == "__main__":
    main()




### example
#python run_compute_cor-disease.py \
    #--input_dir /space/scratch/nairuz-rewiring/data/testing \
    #--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease \
    #--correlation_type pearson \ 
    #--replace_nans True
    #--analysis_type disease \ 
    #--categories AD CTL









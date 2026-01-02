#!~/anaconda3/envs/homl3/bin/python

import sys
import os
from pathlib import Path
import argparse
import scanpy as sc


# Add the path to source the function
sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from computeCor_parallelV1 import compute_gene_correlation
from saveCorDict import save_to_hdf5

#### this is when the script takes in a single file 
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run compute_gene_correlation on a single .h5ad file.")
    parser.add_argument("--correlation_type", type=str, default="pearson", choices=["pearson", "spearman"],
                        help="Type of correlation to compute (default: 'pearson').")
    parser.add_argument("--replace_nans", type=bool, default=True,
                        help="Whether to replace NaNs in the correlation matrix with the median value (default: True).")
    parser.add_argument("--analysis_type", type=str, default=None, choices=["baseline", None],
                        help="Type of analysis. If 'baseline', subset data where diagnosis == 'no'.")
    parser.add_argument("--input_file", type=str, required=True,
                        help="Path to the input .h5ad file.")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Path to the output directory where results will be saved.")
    args = parser.parse_args()


    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check if the input file is valid
    if not input_file.is_file() or input_file.suffix != ".h5ad":
        print(f"Invalid input file: {input_file}")
        sys.exit(1)

    print(f"Processing file: {input_file}")

    try:
        adata = sc.read_h5ad(input_file)
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        sys.exit(1)

    # Define the output file name
    output_file = output_dir / "cor-summary"/ f"{input_file.stem}_corSummary.csv"

    # Skip if already processed
    if output_file.exists():
        print(f"Output already exists for {input_file}. Skipping.")
        return

    # Run the compute_gene_correlation function
    try:
        cor = compute_gene_correlation(
            adata=adata,
            correlation_type=args.correlation_type,
            replace_nans=args.replace_nans,
            min_cells_threshold=20,
            max_workers=6,
            analysis_type=args.analysis_type,
            file_path=input_file
        )
    except Exception as e:
        print(f"Error processing file {input_file}: {e}")
        sys.exit(1)

    # Save the result
    try:
        save_to_hdf5(cor, output_dir, input_file)
        print(f"Saved results to: {output_file}")
    except Exception as e:
        print(f"Error saving results for file {input_file}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()











#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import sys
import argparse
import scanpy as sc
import gc
from pathlib import Path

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV6_disease import compute_gene_correlation

def main():
    parser = argparse.ArgumentParser(description="Compute gene correlation for a single file and category.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to input .h5ad file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument("--category", type=str, required=True, help="Disease category.")
    parser.add_argument("--correlation_type", type=str, default="pearson", choices=["pearson", "spearman"], help="Type of correlation.")
    parser.add_argument("--replace_nans", action='store_true', help="Replace NaNs with median.")

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)

    try:
        adata = sc.read_h5ad(input_file)
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        return

    if args.category not in adata.obs["disease"].unique():
        print(f"Category {args.category} not found in {input_file}. Skipping.")
        return

    try:
        compute_gene_correlation(
            adata=adata,
            correlation_type=args.correlation_type,
            replace_nans=args.replace_nans,
            min_cells_threshold=20,
            category=args.category,
            file_path=input_file,
            output_dir=output_dir
        )
    except Exception as e:
        print(f"Error processing {input_file} for category {args.category}: {e}")

    del adata
    gc.collect()


if __name__ == "__main__":
    main()

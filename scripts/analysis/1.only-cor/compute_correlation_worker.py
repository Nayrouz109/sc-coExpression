#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import sys
import argparse
import scanpy as sc
import gc
from pathlib import Path

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
#from corOnly_parallelV6_disease import compute_gene_correlation
from corOnly_parallelV6_xCellxSubject import compute_gene_correlation
#from corOnly_parallelV6w_disease import compute_gene_correlation
#from corOnly_parallelV6wm_disease import compute_gene_correlation

def main():
    parser = argparse.ArgumentParser(description="Compute gene correlation for a single file and category.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to input .h5ad file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument("--category", type=str, default=None, help="Disease category (or 'ALL' for all samples).")
    parser.add_argument("--correlation_type", type=str, default="pearson",
                        choices=["pearson", "spearman"], help="Type of correlation.")
    parser.add_argument("--replace_nans", action='store_true', help="Replace NaNs with median.")
    parser.add_argument("--xCellxSubject", action='store_true',
                        help="If set, compute correlation across all patients instead of per-patient aggregation.")
    parser.add_argument("--min_cells_threshold", type=int, default=20,
                        help="Minimum number of cells required per patient.")

    args = parser.parse_args()
    input_file = Path(args.input_file)
    output_dir = Path(args.output_dir)

    try:
        adata = sc.read_h5ad(input_file)
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")
        return

    # Handle category logic
    if args.category and args.category not in ["ALL", "None", "none", ""]:
        if args.category not in adata.obs["disease"].unique().tolist():
            print(f"⚠️ Category '{args.category}' not found in {input_file.name}. Skipping.")
            return
        else:
            # Filter to selected category
            adata = adata[adata.obs["disease"] == args.category].copy()
            print(f"✅ Processing {args.category} samples: {adata.n_obs} cells.")
    else:
        # No filtering → use all cells
        print(f"✅ Processing ALL samples together: {adata.n_obs} cells.")

    # Compute correlation
    try:
        compute_gene_correlation(
            adata=adata,
            correlation_type=args.correlation_type,
            replace_nans=args.replace_nans,
            min_cells_threshold=args.min_cells_threshold,
            category=args.category if args.category not in ["ALL", None, "None", "none", ""] else "ALL",
            file_path=input_file,
            output_dir=output_dir,
            xCellxSubject=args.xCellxSubject
        )
    except Exception as e:
        print(f" Error processing {input_file.name} (category={args.category}): {e}")

    # Cleanup
    del adata
    gc.collect()


if __name__ == "__main__":
    main()




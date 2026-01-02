#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

import os
import glob
import sys
import argparse
from pathlib import Path

sys.path.append("/home/nelazzabi/rewiring/scripts")
from functions.technical_controls.process_file import process_file

def main():
    parser = argparse.ArgumentParser(description="Shuffle and run correlation.")
    parser.add_argument("--input_dir", type=str, required=True, help="Path to input .h5ad file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save results.")
    parser.add_argument("--cell_type", type=str, required=True, help="Cell type prefix for filtering files")
    parser.add_argument("--correlation_type", type=str, default="pearson", choices=["pearson", "spearman"], help="Type of correlation.")
    parser.add_argument("--n_shuffles", type=int, required=True, help="Number of shuffling iterations")
    parser.add_argument("--categories", type=str, required=True, help="Comma-separated disease categories (e.g., AD,CTL)")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)  
    output_dir = Path(args.output_dir)

    print(f"Processing files in {input_dir} for cell type: {args.cell_type}")
    
    # Using glob to find all matching files
    for file_path in glob.glob(os.path.join(input_dir, f"{args.cell_type}*.h5ad")):
        process_file(file_path, 
                     set(args.categories.split(",")),  # Ensure categories are passed as a set
                     args.n_shuffles, 
                     output_dir, 
                     args.correlation_type)

if __name__ == "__main__":
    main() 


#python process_files.py 
# --input_dir /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease 
# --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/n-shuffled 
# --cell_type Mic 
# --correlation_type pearson 
# --n_shuffles 2 
# --categories 'AD,CTL'

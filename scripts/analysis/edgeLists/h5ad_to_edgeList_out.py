#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

import os
import anndata as ad
import pandas as pd
import argparse
import sys

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor/edgeList')
from h5ad_to_edgeList import process_h5ad_file

def main():
    parser = argparse.ArgumentParser(description="Process .h5ad files into edge lists.")
    parser.add_argument("input_dir", help="Path to the input directory containing .h5ad files")
    parser.add_argument("output_dir", help="Path to the output directory where results will be saved")
    args = parser.parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process each .h5ad file in the input directory
    for file_name in os.listdir(args.input_dir):
        if file_name.endswith(".h5ad") and file_name.startswith("Mic"):
            input_path = os.path.join(args.input_dir, file_name)
            process_h5ad_file(input_path, args.output_dir)
    
if __name__ == "__main__":
    main()


#python h5ad_to_edgeList_out.py /cosmos/data/project-data/NW-rewiring/coExpr/disease/agg-mtx /cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list
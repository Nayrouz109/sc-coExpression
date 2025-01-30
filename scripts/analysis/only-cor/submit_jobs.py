#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import sys
import os
from pathlib import Path
import argparse
import scanpy as sc
import subprocess

def submit_job(input_file, category, output_dir, correlation_type, replace_nans):
    """
    Submits a SLURM job for computing gene correlation on a given category.
    """
    job_name = f"cor_{Path(input_file).stem}_{category}"
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={output_dir}/logs/{job_name}.out
#SBATCH --error={output_dir}/logs/{job_name}.err
#SBATCH --cpus-per-task=20

source ~/anaconda3/bin/activate homl3

python compute_correlation_worker.py --input_file "{input_file}" --category "{category}" \\
    --output_dir "{output_dir}" --correlation_type "{correlation_type}" \\
    {"--replace_nans" if replace_nans else ""}
"""
    
    # Save SLURM script
    slurm_script_path = f"{output_dir}/slurm_scripts/{job_name}.sh"
    os.makedirs(f"{output_dir}/slurm_scripts", exist_ok=True)
    os.makedirs(f"{output_dir}/logs", exist_ok=True)
    
    with open(slurm_script_path, "w") as f:
        f.write(slurm_script)

    # Submit job
    subprocess.run(["sbatch", slurm_script_path])


def main():
    parser = argparse.ArgumentParser(description="Submit SLURM jobs for compute_gene_correlation.")
    parser.add_argument("--input_dir", type=str, required=True, help="Path to the input directory containing .h5ad files.")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to the output directory for results.")
    parser.add_argument("--correlation_type", type=str, default="pearson", choices=["pearson", "spearman"], help="Type of correlation.")
    parser.add_argument("--replace_nans", action="store_true", help="Replace NaNs in the correlation matrix.")
    parser.add_argument("--categories", type=str, nargs='+', default=['AD', 'CTL'], help="Categories for disease analysis.")

    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for file_path in input_dir.glob("Inh*.h5ad"):
        try:
            adata = sc.read_h5ad(file_path)
            available_categories = set(adata.obs["disease"].unique())
        except Exception as e:
            print(f"Skipping {file_path} due to error: {e}")
            continue

        for category in args.categories:
            if category in available_categories:
                submit_job(str(file_path), category, str(output_dir), args.correlation_type, args.replace_nans)


if __name__ == "__main__":
    main()


### actual directories
# python submit_jobs.py \
    #--input_dir  /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease \
    #--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease \
    #--correlation_type pearson 
    #--replace_nans

### testing directories 
# python submit_jobs.py \
    #--input_dir /space/scratch/nairuz-rewiring/data/testing/test \
    #--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/testing \
    #--correlation_type pearson 
    #--replace_nans
### #SBATCH --time=2:00:00
### #SBATCH --mem=20G
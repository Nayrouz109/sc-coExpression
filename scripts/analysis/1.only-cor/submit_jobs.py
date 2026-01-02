#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import sys
import os
from pathlib import Path
import argparse
import scanpy as sc
import subprocess

def submit_job(input_file, category, output_dir, correlation_type, replace_nans, xcellxsubject=False, min_cells_threshold=None):

    """
    Submits a SLURM job for computing gene correlation on a given category.
    """
    job_name = f"cor_{Path(input_file).stem}_{category}"
    category_arg = f'--category "{category}"' if category not in [None, "ALL"] else ""
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={output_dir}/logs/{job_name}.out
#SBATCH --error={output_dir}/logs/{job_name}.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8


source ~/anaconda3/bin/activate homl3

python /home/nelazzabi/rewiring/scripts/analysis/1.only-cor/compute_correlation_worker.py --input_file "{input_file}" {category_arg} \\
    --output_dir "{output_dir}" --correlation_type "{correlation_type}" \\
    {"--replace_nans" if replace_nans else ""} \\
    {"--xCellxSubject" if xcellxsubject else ""} \\
    {f"--min_cells_threshold {min_cells_threshold}" if min_cells_threshold is not None else ""}


"""
    #SBATCH --cpus-per-task=20
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
    parser.add_argument("--categories", type=str, nargs="*", default=None, help="Categories for disease analysis. If omitted, the correlation will be computed using all samples combined.")
    parser.add_argument("--cell_type", type=str, required=True, help="Cell type to filter .h5ad files (e.g., 'Mic', 'Exc', etc.).")
    parser.add_argument("--xCellxSubject", action="store_true", help="Compute xCellxSubject correlations (filtering patients with too few cells).")
    parser.add_argument("--min_cells_threshold", type=int, default=None, help="Minimum number of cells per patient (used if --xCellxSubject).")

    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    file_pattern = f"{args.cell_type}*.h5ad"
    #for file_path in input_dir.glob("Exc*.h5ad"):

    for file_path in input_dir.glob(file_pattern):
        try:
            adata = sc.read_h5ad(file_path)
            available_categories = set(adata.obs["disease"].unique())
        except Exception as e:
            print(f"Skipping {file_path} due to error: {e}")
            continue

        # 🔹 Case 1: Run for all samples (no category filtering)
        if not args.categories:
            print(f"Submitting job for ALL categories in {file_path.name}")
            submit_job(
                str(file_path),
                "ALL",  # or you could pass None
                str(output_dir),
                args.correlation_type,
                args.replace_nans,
                args.xCellxSubject,
                args.min_cells_threshold
            )
            continue

        # 🔹 Case 2: Run separately per specified category
        for category in args.categories:
            if category in available_categories:
                submit_job(
                    str(file_path),
                    category,
                    str(output_dir),
                    args.correlation_type,
                    args.replace_nans,
                    args.xCellxSubject,
                    args.min_cells_threshold
                )


if __name__ == "__main__":
    main()


### actual directories
# python submit_jobs.py \
    #--input_dir  /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease \
    #--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease \
    #--correlation_type pearson 
    #--replace_nans
    #--cell_type Mic
    #--categories CTL SCZ

### Sep 11 2025 
#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject \
#  --correlation_type pearson \
#  --replace_nans \
#  --cell_type Mic \
#  --categories CTL \
#  --xCellxSubject \
#  --min_cells_threshold 20

### Pans: 1. xCell - Sep 30th, 2025 
#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/cpm \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/0all \
#  --correlation_type pearson \
#  --replace_nans \
#  --cell_type Mic \
#  --categories CTL \ # or leave it out if you want to run it for all samples
#  --min_cells_threshold 20 

### Pans: 2. xCellxSubject 
#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/cpm \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject \
#  --correlation_type pearson \
#  --replace_nans \
#  --cell_type Mic \
#  --categories CTL \
#  --xCellxSubject \
#  --min_cells_threshold 20


















### union genes April 2025 
# python submit_jobs.py --input_dir /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/union-lessStringent 
# --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/0all-union 
# --correlation_type pearson 
# --replace_nans 
# --cell_type Mic  
# --categories CTL AD


### weighted mean 
# python submit_jobs.py \
    #--input_dir  /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease \
    #--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease/wm \         
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

#--input_dir /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/01shuffledDisease/AD-Ctl
#--output_dir /cosmos/data/project-data/NW-rewiring/coExpr/shuffled 

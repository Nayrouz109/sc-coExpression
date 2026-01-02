#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import sys
import os
import glob
from pathlib import Path
import argparse
import scanpy as sc
import subprocess


def submit_job(file_list, file_path, category, output_dir, correlation_type, replace_nans, xType=False, PB=False, min_cells_threshold=None):

    """
    Submits a SLURM job for computing gene correlation on a given category.
    """
    dataset_name = Path(file_path).stem.replace("allCellTypes_", "")
    category_arg = f'--category "{category}"' if category not in [None, "ALL"] else ""
    job_name = f"cor_{dataset_name}_{category}"
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={output_dir}/logs/{job_name}.out
#SBATCH --error={output_dir}/logs/{job_name}.err
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8


source ~/anaconda3/bin/activate homl3

python /home/nelazzabi/rewiring/scripts/analysis/1.only-cor/xSubjxType/worker.py \
    --file_list { " ".join(file_list) } \
    --file_path "{file_path}" \
    --output_dir "{output_dir}" \
    {category_arg} \
    --correlation_type "{correlation_type}" \
    {"--replace_nans" if replace_nans else ""} \
    {"--xType" if xType else ""} \
    {"--PB" if PB else ""} \
    {f"--min_cells_threshold {min_cells_threshold}" if min_cells_threshold is not None else ""}



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
    #parser.add_argument("--categories", type=str, nargs='+', default=['AD', 'CTL'], help="Categories for disease analysis.")
    parser.add_argument("--categories", type=str, nargs="*", default=None, help="Categories for disease analysis. If omitted, the correlation will be computed using all samples combined.")
    parser.add_argument("--xType", action='store_true', help="If set, compute correlation across all patients and cell types.")
    parser.add_argument("--min_cells_threshold", type=int, default=None, help="Minimum number of cells per patient (used if --xCellxSubject).")
    parser.add_argument("--PB", action='store_true', help="If set, generate PB matrices.") 
    
    args = parser.parse_args()
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    dataset = [
    'ling-2024', 'ROSMAP-Mathys-2024', 'MSBB-2024', 'RADC-2024', 'HBCC-2024', 'SZBDMulti-2024',
    'MSSM-2024', 'Rexach-2024', 'ROSMAP-erosion-2023', 'Kamath-NPH-2023',
    'Hoffman-2023', 'ROSMAP-Mathys-2019', 'ROSMAP-Mathys-2023', 'ROSMAP-TREM2-2020',
    'ROSMAP-DeJager-2024', 'ROSMAP-DeJager-2023', 'SEA-AD-MTG-2024', 'SEA-AD-DLPFC-2024', 
    'Lau-2020', 'Lim-2022', 'Nagy-2020', 'Pineda-2021', 'Velmeshev-2019'
    ]

    for ds in dataset:
        pattern = f"*{ds}.h5ad"
        matching_files = glob.glob(os.path.join(input_dir, pattern))
        
        if not matching_files:
            print(f"Dataset: {ds} — no matching files found.")
            continue

        file_path = os.path.join(input_dir, f"allCellTypes_{ds}.h5ad")

        # 🔹 CASE 1: No categories specified or explicitly "ALL"
        if not args.categories or (len(args.categories) == 1 and args.categories[0].upper() == "ALL"):
            print(f"Submitting job for {ds} — ALL samples combined.")
            submit_job(
                matching_files,
                str(file_path),
                "ALL",                     # category flag for clarity
                str(output_dir),
                args.correlation_type,
                args.replace_nans,
                args.xType,
                args.PB,
                args.min_cells_threshold
            )
            continue

        # 🔹 CASE 2: Run separately for specified categories
        for category in args.categories:
            print(f"Submitting job for {ds} — Category: {category}")
            submit_job(
                matching_files,
                str(file_path),
                category,
                str(output_dir),
                args.correlation_type,
                args.replace_nans,
                args.xType,
                args.PB,
                args.min_cells_threshold
            )



if __name__ == "__main__":
    main()


##############################################################################################################
############################################ xSubjectxType ################################################### 
### Sep 24th 2025 
#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/raw \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType \
#  --correlation_type pearson \
#  --replace_nans \
#  --categories CTL \
#  --xType \
#  --PB \
#  --min_cells_threshold 20


### Sep 29thth 2025 
#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/raw \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType \
#  --correlation_type pearson \
#  --replace_nans \
#  --categories CTL \
#  --xType \
#  --PB \
#  --min_cells_threshold 20

##############################################################################################################
############################################ xSubjectx ####################################################### 
### Sep 30th 2025 
#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/raw \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/xSubject \
#  --correlation_type pearson \
#  --replace_nans \
#  --categories CTL \
#  --PB \
#  --min_cells_threshold 20

#python submit_jobs.py \
#  --input_dir /cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/raw \
#  --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/xSubject \
#  --correlation_type pearson \
#  --replace_nans \
#  --categories CTL \
#  --PB \
#  --min_cells_threshold 20
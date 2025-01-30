#!/usr/bin/bash
#SBATCH --job-name=scCor_sngl            # Job name
#SBATCH --output=logs/%x_%j.out          # Standard output log
#SBATCH --error=logs/%x_%j.err           # Error log
#SBATCH --ntasks=1                       # Single task
#SBATCH --cpus-per-task=10               # Number of CPU cores per task


source ~/anaconda3/bin/activate homl3


# Run the Python script for a single file
python /home/nelazzabi/rewiring/scripts/analysis/run_compute_cor.py \
    --correlation_type pearson \
    --replace_nans True \
    --analysis_type baseline \
    --input_file /space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag/Mic_ROSMAP-DeJager-2024.h5ad \
    --output_dir /space/scratch/nairuz-rewiring/coExpr




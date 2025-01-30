#!/usr/bin/bash
#SBATCH --job-name=scCor_single          # Job name
#SBATCH --output=logs/%x_%j.out          # Standard output log
#SBATCH --error=logs/%x_%j.err           # Error log
#SBATCH --ntasks=1                       # Single task
#SBATCH --cpus-per-task=7                # 6 cores for each file
#SBATCH --mail-user=nayrouz@student.ubc.ca # Email for notifications
#SBATCH --mail-type=ALL                  # Notifications on job events

# Activate your conda environment
source ~/anaconda3/bin/activate homl3

# Input file (set as a variable or pass as an argument to this script)
FILE_TO_PROCESS=$1

if [[ -z "$FILE_TO_PROCESS" ]]; then
    echo "Error: No input file specified."
    exit 1
fi

# Run the Python script for a single file
python /home/nelazzabi/rewiring/scripts/analysis/run_compute_cor.py \
    --correlation_type pearson \
    --replace_nans True \
    --analysis_type baseline \
    --input_file "$FILE_TO_PROCESS" \
    --output_dir /space/scratch/nairuz-rewiring/coExpr

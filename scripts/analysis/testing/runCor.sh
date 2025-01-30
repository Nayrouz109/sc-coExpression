#!/usr/bin/bash
#SBATCH --job-name=scCor_batch             # Job name
#SBATCH --output=logs/%x_%A_%a.out         # Standard output and error log
#SBATCH --error=logs/%x_%A_%a.err          # Error log file
#SBATCH --ntasks=1                         # Single task
#SBATCH --cpus-per-task=18                 # 6 cores per file x 3 files = 18 cores
#SBATCH --array=0-35                       # Adjust range based on batch count
#SBATCH --mail-user=nayrouz@student.ubc.ca # Email for notifications
#SBATCH --mail-type=ALL                    # Notifications on job events

# Activate your conda environment
source ~/anaconda3/bin/activate homl3

# Get the batch file based on array task ID
BATCH_FILE=$(ls batch_* | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

# Loop through each file in the batch and process in parallel
while read -r FILE; do
  python /home/nelazzabi/rewiring/scripts/analysis/run_compute_cor.py \
    --correlation_type pearson \
    --replace_nans True \
    --analysis_type baseline \
    --input_dir "$(dirname "$FILE")" \
    --output_dir /space/scratch/nairuz-rewiring/coExpr &
done < "$BATCH_FILE"

# Wait for all background tasks to complete
wait

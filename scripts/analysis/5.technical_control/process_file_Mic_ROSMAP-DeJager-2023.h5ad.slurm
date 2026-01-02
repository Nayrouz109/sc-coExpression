#!/bin/bash
#SBATCH --job-name=cor_Mic_ROSMAP-DeJager-2023
#SBATCH --output=/cosmos/data/project-data/NW-rewiring/coExpr/n-shuffled/logs/cor_Mic_ROSMAP-DeJager-2023.out
#SBATCH --error=/cosmos/data/project-data/NW-rewiring/coExpr/n-shuffled/logs/cor_Mic_ROSMAP-DeJager-2023.err
#SBATCH --cpus-per-task=20


# Run the Python script for this file
python /home/nelazzabi/rewiring/scripts/functions/technical_controls/process_file.py --file_path /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/00testing/Mic_ROSMAP-DeJager-2023.h5ad  --categories ('AD', 'CTL') --n_shuffles 2 --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/n-shuffled --correlation_type pearson

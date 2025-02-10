#!~/anaconda3/bin/Rscript

input_dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list"
output_dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/disease/top200"
slurm_dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/disease/slurm-scripts"


if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(slurm_dir)) dir.create(slurm_dir, recursive = TRUE)


input_files <- list.files(path = input_dir, pattern = "Mic_.*_corEdgeLists.csv", full.names = TRUE)

# create SLURM script for each file 
for (file in input_files) {
  base_name <- tools::file_path_sans_ext(basename(file))
  slurm_script <- file.path(slurm_dir, paste0(base_name, ".slurm"))
  
  script_content <- sprintf("#!/bin/bash
#SBATCH --job-name=%s
#SBATCH --output=%s/%s.out
#SBATCH --error=%s/%s.err
#SBATCH --cpus-per-task=2

Rscript /home/nelazzabi/rewiring/scripts/functions/topCoExpr/process_top200.r %s %s
", base_name, slurm_dir, base_name, slurm_dir, base_name, file, output_dir)
  
  writeLines(script_content, slurm_script)
  system(paste("sbatch", slurm_script))
}

print("SLURM scripts generated and jobs submitted!")

#/space/opt/bin/Rscript generate_slurm_scripts.r
#cd /cosmos/data/project-data/NW-rewiring/coExpr/disease/slurm-scripts
#for script in *.slurm; do sbatch "$script"; done
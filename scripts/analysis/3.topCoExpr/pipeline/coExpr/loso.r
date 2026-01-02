
#!~/anaconda3/envs/biof-r/bin/R

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: script.R <input_dir> <cell_type> <output_dir>")
}

input_dir <- args[1]
cell_type <- args[2]
output_dir <- args[3]

input_dir <- args[1]
cell_type <- args[2]
output_dir <- args[3]


# Get input files
input_files <- list.files(path = input_dir, pattern = paste0(cell_type, ".*_edgeList.csv"), full.names = TRUE)
cat("Found", length(input_files), "files to process.\n")
TFs <- readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS")
TFs <- lapply(TFs, function(df) df %>% mutate(geneB = Symbol))

# Source necessary functions
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")

# Extract top 200 partners per TF
cat("Extracting top 200 partners per TF...\n")
tf_top200_list <- setNames(
  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = 200)),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)

# Generate TF reproducibility report
cat("Generating TF reproducibility report...\n")
tf_reprod <- calculate_tf_target_reproducibility_loso(tf_top200_list)
output_file3 <- file.path(output_dir, paste0(cell_type, "_tf_reprod_tbl_loso.rds"))
cat("Saving TF reproducibility report to:", output_file3, "\n")
saveRDS(tf_reprod, output_file3)


# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list/all-matrix
# Mic
# /home/nelazzabi/rewiring/data/coExpr/AD-Ctl 




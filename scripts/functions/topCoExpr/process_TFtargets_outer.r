#!~/anaconda3/envs/biof-r/bin/R

library(data.table)
library(dplyr)
library(argparse)

parser <- ArgumentParser(description = "Process to filter top coExpr")
parser$add_argument("input_dir", help = "path to edge list files")
parser$add_argument("output_dir", help = "path to output dir")

args <- parser$parse_args()

# Ensure output directory exists
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
}


source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")
TFs <- readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS") 
input_files <- list.files(path = args$input_dir , pattern = "Mic_.*_edgeList.csv", full.names = TRUE)



tf_top200_list <- setNames(
  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = 200)),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)

tf_allTrg_list <- setNames(
  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = NULL)),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)



output_file1 <- file.path(args$output_dir, "tf_top200_list.rds")
output_file2 <- file.path(args$output_dir, "tf_allTrg_list.rds")

saveRDS(tf_top200_list, output_file1)
saveRDS(tf_allTrg_list, output_file2)





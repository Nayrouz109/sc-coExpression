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
condition <- args[4]

cat("\nProcessing cell type:", cell_type, "\n")

# Get input files
input_files <- list.files(path = input_dir, pattern = paste0(cell_type, ".*_" , condition, ".*_edgeList.csv"), full.names = TRUE)
cat("Found", length(input_files), "files to process.\n")

# Load additional data
TFs <- readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS")
TFs <- lapply(TFs, function(df) df %>% mutate(geneB = Symbol))
tf_by_gene_binding <- readRDS("/home/nelazzabi/rewiring/data/modalities/gene_by_tf_binding_ranked_scaled.rds")
lit_targets = readRDS("/home/nelazzabi/rewiring/data/modalities/lit_curated_targets.rds")


# Source necessary functions
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")

# 0. Extract all partners per TF
cat("Extracting all partners per TF...\n")
tf_allTrg_list <- setNames(
  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = NULL)),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)

for (i in names(tf_allTrg_list)) {
  tf_allTrg_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
}


#output_tf_allTrg <- file.path(output_dir, paste0(cell_type, "_tf_allTrg_list.rds"))
#cat("Saving all TF target list to:", output_tf_allTrg, "\n")
#saveRDS(tf_allTrg_list, output_tf_allTrg) 


# aggregate TF-gene coexpression profiles and add: 
# binding score 
# literature evidence 

summary_tf_allTargets <- summarize_TF_target_profiles(tf_allTrg_list, condition, tf_by_gene_binding, lit_targets)

output <- file.path(output_dir, paste0(cell_type, "_tf_allTrg_avgRnkCoexpr_binding_lit.rds"))
cat("Saving TF target binding and lit:", output, "\n")
saveRDS(summary_tf_allTargets, output) 







# 1. calculate average ranking for each TF-gene pair 
#cat("calculate average ranking for each TF-gene pair...\n")
#summary_tf_allTargets <- average_tf_target_ranks(tf_allTrg_list, condition, tf_by_gene_binding)

#output_tf_allTrg <- file.path(output_dir, paste0(cell_type, "_tf_allTrg_avgRnkCoexpr_binding.rds"))
#cat("Saving TF target binding:", output_tf_allTrg, "\n")
#saveRDS(summary_tf_allTargets, output_tf_allTrg) 

# 2. generate TF-gene coexpression aggregate profiles 
#cat("generate TF-gene coexpression aggregate profiles...\n")
#summary_tf_allTargets <- tf_gene_aggCoExpr_litCur(tf_allTrg_list, condition, tf_by_gene_binding)

#output_tf_allTrg <- file.path(output_dir, paste0(cell_type, "_tf_allTrg_avgRnkCoexpr_litCur.rds"))
#cat("Saving TF target lit curation:", output_tf_allTrg, "\n")
#saveRDS(summary_tf_allTargets, output_tf_allTrg) 


# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/0all/edge-list/all-matrix
# Mic
# /home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell
# CTL 



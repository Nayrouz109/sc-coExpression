
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
output_dir <- args[2]
condition <- args[3]


# Load additional data
TFs <- readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS")
TFs <- lapply(TFs, function(df) df %>% mutate(geneB = Symbol))
tf_by_gene_binding <- readRDS("/home/nelazzabi/rewiring/data/modalities/gene_by_tf_binding_ranked_scaled.rds")
lit_targets = readRDS("/home/nelazzabi/rewiring/data/modalities/lit_curated_targets.rds")
tf_rpr = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell/intgr/TF_rpr_inegrated_cellTypes.rds") %>% setDT()


# Source necessary functions
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")


# 1. Read top 200 partners per TF
cat("reading top 200 partners per TF...\n")

dir_path <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell"
rds_files <- list.files(path = dir_path, pattern = "tf_top200_list\\.rds$", full.names = TRUE)

tf_top200_list <- setNames(
  lapply(rds_files, readRDS),
  gsub("_tf_top200_list\\.rds$", "", basename(rds_files))
)




# 2. Loop through all cell types in tf_top200_list
for (cell_type in names(tf_top200_list)) {
  cat("Processing cell type:", cell_type, "\n")
  
  # Apply summarization function
  result <- summarize_TF_target_profilesV1(
    tf_top200_list = tf_top200_list,
    condition = condition,
    cell_type = cell_type,
    tf_by_gene_binding = tf_by_gene_binding,
    lit_targets = lit_targets,
    tf_rpr = tf_rpr
  )
  
  # Construct output file path
  output_file <- file.path(output_dir, paste0(cell_type, "_tf_200Trg_avgRnkCoexpr_binding_lit_TFrpr.rds"))
  
  # Save result and log progress
  cat("Saving TF top 200 target + binding + lit + rpr score for cell type [", cell_type, "]:", output_file, "\n")
  saveRDS(result, output_file)
}


# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/0all/edge-list/all-matrix
# /home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell
# CTL 




### -----------------------------------------------------------------------------------------------------------------------------------
### this is to summarize the freq per gene pair and generate 1 summary dataframe 
### --------------------------------------------------------------------------------------------------------------------------------------
# Step 1: List relevant files
rds_files <- list.files(
  path = output_dir,
  pattern = "_tf_200Trg_avgRnkCoexpr_binding_lit_TFrpr\\.rds$",
  full.names = TRUE
)

# Step 2: Read all into a list and combine into a data.table
print("reading the files i just saved") 
tf_target_list <- lapply(rds_files, readRDS)
tf_target_dt <- rbindlist(tf_target_list)

# Step 3: Create unique key for each geneA-geneB pair
tf_target_dt[, pair_key := paste(geneA, geneB, sep = "_")]

# Step 4: Count frequency of each pair
pair_counts <- tf_target_dt[, .(pair_freq = .N), by = pair_key]

# Step 5: Merge counts back into original data
print("merging the counts with the data summary") 
tf_target_dt <- merge(tf_target_dt, pair_counts, by = "pair_key")

# Step 6: adding target avg expression 
allGenes_xpr = fread("/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/allGenes_exp_summary_logCPM.csv") %>% setDT()
allGenes_avg <- allGenes_xpr[, .(target_avg_expr_log_cpm = mean(avg_expr_log_cpm, na.rm = TRUE)), 
                             by = .(gene, cell_type)]

# step 2: join with tf_target_dt (geneB <-> gene, and cell_type match)
tf_target_dt <- merge(
  tf_target_dt, 
  allGenes_avg, 
  by.x = c("geneB", "cell_type"), 
  by.y = c("gene", "cell_type"), 
  all.x = TRUE
)


# Step 7: Optional: reorder columns (put pair_freq last or wherever you prefer)
setcolorder(tf_target_dt, c(
  "geneA", "geneB", "avg_rank_norm_corr", "n_datasets", "condition", "cell_type",
  "binding_score", "TF_target_predicted", "tf_rpr_score", "pair_freq"
))

cat("Saving TF top 200 target + binding + lit + rpr score + freq", "\n")
output_file <- file.path(output_dir, paste0("intg_TF_200Trg_avgRnkCoexpr_binding_lit_TFrpr.rds"))
saveRDS(tf_target_dt, output_file)





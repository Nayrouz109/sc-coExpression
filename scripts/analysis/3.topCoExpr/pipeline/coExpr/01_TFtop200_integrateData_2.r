#!~/anaconda3/envs/biof-r/bin/R




#### this is only the second bottom half of 01_TFtop200_integrateData.r script
# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
})


output_dir = "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell"
# Step 1: List relevant files
rds_files <- list.files(
  path = "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell",
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
pair_counts <- tf_target_dt[, .(
  pair_freq = .N, 
  cell_types_present = paste(unique(cell_type), collapse = ",")
), by = pair_key]


# Step 5: Merge counts back into original data
print("merging the counts with the data summary") 
tf_target_dt <- merge(tf_target_dt, pair_counts, by = "pair_key")

# Step 6: adding target avg expression 
allGenes_xpr = fread("/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/allGenes_exp_summary_logCPM.csv") %>% setDT()
allGenes_avg <- allGenes_xpr[, .(target_avg_expr_log_cpm = mean(avg_expr_log_cpm, na.rm = TRUE)), 
                             by = .(gene, cell_type)]

tf_target_dt <- merge(
  tf_target_dt, 
  allGenes_avg, 
  by.x = c("geneB", "cell_type"), 
  by.y = c("gene", "cell_type"), 
  all.x = TRUE
)

# Step 7: adding target avg expression 
TF_chipseq_nDataset = readRDS("/home/nelazzabi/rewiring/data/modalities/TF_chipseq_ndatasets.rds")
tf_target_dt <- merge(
  tf_target_dt,
  TF_chipseq_nDataset,
  by = "geneA",
  all.x = TRUE
)

# Step 8: Optional: reorder columns (put pair_freq last or wherever you prefer)
setcolorder(tf_target_dt, c(
  "geneA", "geneB", "avg_rank_norm_corr", "n_datasets", "condition", "cell_type",
  "binding_score", "TF_target_predicted", "tf_rpr_score", "pair_freq"
))

cat("Saving TF top 200 target + binding + lit + rpr score + freq", "\n")
output_file <- file.path(output_dir, paste0("intg_TF_200Trg_avgRnkCoexpr_binding_lit_TFrpr.rds"))
saveRDS(tf_target_dt, output_file)





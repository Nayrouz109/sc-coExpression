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

cat("\nProcessing cell type:", cell_type, "\n")

# Get input files
input_files <- list.files(path = input_dir, pattern = paste0(cell_type, ".*_edgeList.csv"), full.names = TRUE)
cat("Found", length(input_files), "files to process.\n")

# Load additional data
TFs <- readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS")
TFs <- lapply(TFs, function(df) df %>% mutate(geneB = Symbol))
lit_TFs <- fread("/home/nelazzabi/rewiring/data/TFs/Curated_targets_all_Sept2023.tsv")
eRegulon <- fread("/home/nelazzabi/rewiring/data/TFs/eRegulon_metadata_filtered.csv")
microglia_TFs_chatGPT <- c("SPI1", "IRF8", "SALL1", "MEF2C", "MafB")
microglia_TFs_paper <- colnames(fread("/home/nelazzabi/rewiring/data/TFs/cellType-TFs/mic/Supplemental_table_3.csv"))[-1]
mic_TF <- unique(c(microglia_TFs_chatGPT, microglia_TFs_paper))
mic_neuroinflammation <- fread("/home/nelazzabi/rewiring/data/AD-genes/AD_genes_mic_neuroinflammation.csv")
brain_TFs <- fread("/home/nelazzabi/rewiring/data/TFs/brain_tfs.tsv", select = "TF_Symbol")
brain_TF_list <- brain_TFs$TF_Symbol
SEA_AD_micGenes <- fread("/home/nelazzabi/rewiring/data/AD-genes/SEAAD_mic_genes.csv")
ribo = fread("/home/nelazzabi/rewiring/data/genes/ribosomal/group-1054.csv", header = TRUE) %>% filter(Group %in% c("L ribosomal proteins", "S ribosomal proteins")) 
ribo_list <- split(ribo, ribo$`Approved symbol`)  
names(ribo_list) <- make.unique(names(ribo_list))
cat("done reading additional files...\n")

# Source necessary functions
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")


# 0.1 ribosomal reproducability 
cat("Extracting top 89 partners per ribo...\n")
ribo_top_list <- setNames(
  lapply(input_files, function(file) process_file(file, ribo_list, top_n_targets = length(ribo_list))),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)

for (i in names(ribo_top_list)) {
  ribo_top_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
}
cat("Generating ribo reproducibility report...\n")
ribo_reprod <- calculate_tf_target_reproducibility(ribo_top_list)
output_file_ribo <- file.path(output_dir, paste0(cell_type, "_t" ,length(ribo_list), "_rpr_rb.rds"))
cat("Saving ribo reproducibility report to:", output_file_ribo, "\n")
saveRDS(ribo_reprod, output_file_ribo) 

# 1. Extract top 200 partners per TF
cat("Extracting top 200 partners per TF...\n")
tf_top200_list <- setNames(
  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = 200)),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)
for (i in names(tf_top200_list)) {
  tf_top200_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
}

output_tfTop200 <- file.path(output_dir, paste0(cell_type, "_tf_top200_list.rds"))
cat("Saving top 200 TF target list to:", output_tfTop200, "\n")
saveRDS(tf_top200_list, output_tfTop200)

# 0.2 generating a TF null 
cat("generating TF null...\n")
TF_null <- compute_null_overlap(tf_top200_list, 1000) 
#mutate(geneA = paste0("gene_", iteration), meanTop200Ctl = avg_overlap)
output_file_null <- file.path(output_dir, paste0(cell_type, "_t200_rpr_tn.rds"))
cat("Saving TF null reproducibility report to:", output_file_null, "\n")
saveRDS(TF_null, output_file_null)

# 2. Generate TF reproducibility report
cat("Generating TF reproducibility report...\n")
tf_reprod <- calculate_tf_target_reproducibility(tf_top200_list)
tf_reprod$mic_TF <- tf_reprod$geneA %in% mic_TF
tf_reprod$brain_TF = tf_reprod$geneA %in% brain_TF_list
output_tf_reprod <- file.path(output_dir, paste0(cell_type, "_t200_rpr_tf.rds"))
cat("Saving TF reproducibility report to:", output_tf_reprod, "\n")
saveRDS(tf_reprod, output_tf_reprod)
cat("Processing completed successfully.\n")

# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/0all/edge-list/all-matrix
# Mic
# /home/nelazzabi/rewiring/data/coExpr/0all-human
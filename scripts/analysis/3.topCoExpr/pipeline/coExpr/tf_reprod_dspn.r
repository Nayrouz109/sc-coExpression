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

# create dsp list 
dsp_files   <- list.files(path = input_dir, pattern = paste0(cell_type, ".*dsp.*_edgeList.csv"), full.names = TRUE)
non_dsp_files <- setdiff(input_files, dsp_files)
#studies <- unique(sapply(strsplit(basename(input_files), "_"), function(x) x[2]))
dsp_iterations <- unique(gsub(".*_dsp([0-9]+)_.*", "\\1", dsp_files))

input_files_list <- vector("list", length = length(dsp_iterations))
names(input_files_list) <- paste0("dsp", 1:length(dsp_iterations))

for (i in 1:length(dsp_iterations)) {
  # 1. Select files for the current iteration
  dsp_subset <- grep(paste0("_dsp", i, "_"), dsp_files, value = TRUE)
  # 2. Combine with non-dsp files (present in all iterations) 
  input_files_list[[i]] <- c(dsp_subset, non_dsp_files)
}


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
cat("done reading additional files...\n")

# Source necessary functions
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")

# 1. Extract top 200 partners per TF
cat("Extracting top 200 partners per TF...\n")
tf_top200_list_all <- lapply(input_files_list, function(input_files) {
  setNames(
    lapply(input_files, function(file) process_file(file, TFs, top_n_targets = 200)),
    gsub("_edgeList\\.csv$", "", basename(input_files))
  )
})

tf_top200_list_all <- lapply(tf_top200_list_all, function(tf_top200_list) {
  for (i in names(tf_top200_list)) {
    tf_top200_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
  }
  return(tf_top200_list)
})


# 2. Generate TF reproducibility report
cat("Generating TF reproducibility report...\n")
tf_reprod_all <- lapply(tf_top200_list_all, function(tf_top200_list) {
  tf_reprod <- calculate_tf_target_reproducibility(tf_top200_list)
  tf_reprod$mic_TF <- tf_reprod$geneA %in% mic_TF
  tf_reprod$brain_TF <- tf_reprod$geneA %in% brain_TF_list
  return(tf_reprod)
})

names(tf_reprod_all) <- names(tf_top200_list_all)
output_file3 <- file.path(output_dir, paste0(cell_type, "_tf_reprod_tbl.rds"))
cat("Saving TF reproducibility report to:", output_file3, "\n")
saveRDS(tf_reprod_all, output_file3)


# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/downsample/AD-Ctl-dsp-n/edge-list/all-matrix
# Mic
# /home/nelazzabi/rewiring/data/coExpr/downsampling/AD-Ctl-dsp-n

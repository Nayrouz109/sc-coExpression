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
#source("/home/nelazzabi/rewiring/scripts/functions/mm_cons/cons_analysis.r")

# 0.1 ribosomal reproducability 
#cat("Extracting top 89 partners per ribo...\n")
#ribo_top_list <- setNames(
#  lapply(input_files, function(file) process_file(file, ribo_list, top_n_targets = length(ribo_list))),
#  gsub("_edgeList\\.csv$", "", basename(input_files))
#)
#for (i in names(ribo_top_list)) {
#  ribo_top_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
#}
#cat("Generating ribo reproducibility report...\n")
#ribo_reprod <- calculate_tf_target_reproducibility(ribo_top_list)
#output_file_ribo <- file.path(output_dir, paste0(cell_type, "_ribo_reprod_tbl.rds"))
#cat("Saving TF reproducibility report to:", output_file_ribo, "\n")
#saveRDS(ribo_reprod, output_file_ribo) 



cat("Extracting all partners per ribo gene...\n")
ribo_allTrg_list <- setNames(
  lapply(input_files, function(file) process_file(file, TFs = ribo_list, top_n_targets = NULL)),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)

for (i in names(ribo_allTrg_list)) {
  ribo_allTrg_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
}

output_file_ribo <- file.path(output_dir, paste0(cell_type, "_ribo_allGenes_allTrg_list.rds"))
cat("Saving all ribo target list to:", output_file_ribo, "\n")
saveRDS(ribo_allTrg_list, output_file_ribo)

# 1. Extract top 200 partners per TF
#cat("Extracting top 200 partners per TF...\n")
#tf_top200_list <- setNames(
#  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = 200)),
#  gsub("_edgeList\\.csv$", "", basename(input_files))
#)
#for (i in names(tf_top200_list)) {
#  tf_top200_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
#}

#output_file1 <- file.path(output_dir, paste0(cell_type, "_tf_top200_list.rds"))
#cat("Saving top 200 TF target list to:", output_file1, "\n")
#saveRDS(tf_top200_list, output_file1)

# 2. Extract all partners per TF
#cat("Extracting all partners per TF...\n")
#tf_allTrg_list <- setNames(
#  lapply(input_files, function(file) process_file(file, TFs, top_n_targets = NULL)),
#  gsub("_edgeList\\.csv$", "", basename(input_files))
#)

#for (i in names(tf_allTrg_list)) {
#  tf_allTrg_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
#}

#output_file2 <- file.path(output_dir, paste0(cell_type, "_tf_allTrg_list.rds"))
#cat("Saving all TF target list to:", output_file2, "\n")
#saveRDS(tf_allTrg_list, output_file2)

# 3. Generate TF reproducibility report
#cat("Generating TF reproducibility report...\n")
#tf_reprod <- calculate_tf_target_reproducibility(tf_top200_list)
#tf_reprod$mic_TF <- tf_reprod$geneA %in% mic_TF
#tf_reprod$brain_TF = tf_reprod$geneA %in% brain_TF_list
#output_file3 <- file.path(output_dir, paste0(cell_type, "_tf_reprod_tbl.rds"))
#cat("Saving TF reproducibility report to:", output_file3, "\n")
#saveRDS(tf_reprod, output_file3)

# 3.1 generating a TF null 
#cat("generating TF null...\n")
#TF_null <- compute_null_overlap(tf_top200_list, 1000) 
#output_file_null <- file.path(output_dir, paste0(cell_type, "_tf_null_tbl.rds"))
#cat("Saving TF reproducibility report to:", output_file_null, "\n")
#saveRDS(TF_null, output_file_null)



# Generate TF-target summary
#cat("Generating TF-target summary...\n")
#freqTable <- summarize_gene_pairs(tf_top200_list)
#tf_top200_cor_intg_list <- list()


#valid_tfs <- unique(freqTable$geneA)
#cormat_FZ <- readRDS("/home/nelazzabi/rewiring/data/TFs/aggregate_cormat_FZ_hg.RDS")[["Agg_mat"]]
#cormat_allRank <- readRDS("/home/nelazzabi/rewiring/data/TFs/aggregate_cormat_allrank_filter_lenient_hg.RDS")[["Agg_mat"]]

#for (tf in valid_tfs) {
#  cat("processing", tf, "\n")

  
#  if (tf %in% rownames(cormat_allRank) & tf %in% rownames(cormat_FZ)) {
#    tf_targets_alex_fishZ <- rank_gene_partners(cormat_FZ, tf)
#    tf_targets_alex_allRank <- rank_gene_partners(cormat_allRank, tf)
      
#    tf_top200_cor_intg_list[[tf]] <- integrate_gene_ranking(
#      freqTable, 
#      tf_targets_alex_fishZ, 
#      tf_targets_alex_allRank, 
#      TFs, 
#      tf, 
#      lit_TFs, 
#      NULL, 
#      mic_neuroinflammation, 
#      SEA_AD_micGenes, 
#      eRegulon
#    )

#  } else {
#    tf_top200_cor_intg_list[[tf]] <- NA
#  }
  
#  cat("done processing", tf, "\n") 
#}


#output_file4 <- file.path(output_dir, paste0(cell_type, "_tf_top200_cor_intg_list.rds"))
#saveRDS(tf_top200_cor_intg_list, output_file4)
#cat("Processing completed successfully.\n")



# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list/all-matrix
# Mic
# /home/nelazzabi/rewiring/data/coExpr/AD-Ctl 
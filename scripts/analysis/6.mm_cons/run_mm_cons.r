#!~/anaconda3/envs/biof-r/bin/R

# Load required libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: script.R <input_dir> <cell_type> <condition> <output_dir> <mouse_input> <mm_agg>")
}

input_dir <- args[1]
cell_type <- args[2]
condition <- args[3]
output_dir <- args[4]
mouse_input <- args[5] 
mm_agg <- args[6] 

cat("\nProcessing cell type:", cell_type, "\n")


# source functions 
source("/home/nelazzabi/rewiring/scripts/functions/mm_cons/cons_analysis.r")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")

# DATA ==================================================================================================
# =======================================================================================================
# 1. Get input files - human 
input_files <- list.files(path = input_dir, pattern = paste0(cell_type, ".*_edgeList.csv"), full.names = TRUE)
cat("Found", length(input_files), "files to process.\n")

# 2. Get input data - mouse 
# Mouse and human all rank (using the filtering Alex presented)
cat("reading mouse cor matrix .\n")
alex_mouse = readRDS(mouse_input)

genes_map = fread("/home/nelazzabi/rewiring/data/genes/mouse_human/hg_mm_1to1_ortho_genes_DIOPT_V9.tsv")
mouse_to_human <- setNames(genes_map$Symbol_hg, genes_map$Symbol_mm)

flt_mm <- convert_mouse_to_human(alex_mouse$Agg_mat, mouse_to_human)
rm(alex_mouse)

# 3. match mouse and human genes 
cat("matching mm and hg matrices .\n")
matched_data <- filter_edge_lists_and_matrix(input_files, flt_mm)
matched_data_hg <- matched_data$filtered_edges
matched_data_mm <- matched_data$filtered_matrix

cat("converting mm into edge list .\n")
matched_data_mm = generate_edge_list(matched_data_mm, top_n_targets = 200)
rm(matched_data)
rm(flt_mm) 

# 4. Load additional data
TFs <- readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS")
TFs <- lapply(TFs, function(df) df %>% mutate(geneB = Symbol))

ribo = fread("/home/nelazzabi/rewiring/data/genes/ribosomal/group-1054.csv", header = TRUE) %>% filter(Group %in% c("L ribosomal proteins", "S ribosomal proteins")) 
ribo_list <- split(ribo, ribo$`Approved symbol`)  
names(ribo_list) <- make.unique(names(ribo_list))
cat("done reading additional files...\n")





# 1. Extract all partners per gene ==================================================================================================
# =================================================================================================================================
cat("Extracting all partners per gene...\n")

gn_allTrg_list <- setNames(
  lapply(matched_data_hg, function(file) process_file(file, TFs = NULL, top_n_targets = NULL)),
  gsub("_edgeList\\.csv$", "", names(matched_data_hg))
)

for (i in names(gn_allTrg_list)) {
  gn_allTrg_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
}

output_file1 <- file.path(output_dir, paste0(cell_type, "_", mm_agg, "_allGenes_allTrg_list.rds"))
#cat("Saving all gene targets list to:", output_file1, "\n")
#saveRDS(gn_allTrg_list, output_file1)

# 2. Extract top 200 partners per gene ==================================================================================================
# =================================================================================================================================
cat("Extracting top 200 partners per gene...\n")
gn_top200_list <- setNames(
  lapply(matched_data_hg, function(file) process_file(file, TFs = NULL, top_n_targets = 200)),
  gsub("_edgeList\\.csv$", "", names(matched_data_hg))
)
for (i in names(gn_top200_list)) {
  gn_top200_list[[i]]$condition <- ifelse(grepl("_AD", i), "AD", "control")
}

output_file2 <- file.path(output_dir, paste0(cell_type, "_", mm_agg, "_gn_top200_list.rds"))
#cat("Saving top 200 gn target list to:", output_file2, "\n")
#saveRDS(gn_top200_list, output_file2)

rm(matched_data_hg)
# 3. gene-target freq table =================================================================================================== 
# =================================================================================================================================
cat("Generating gene-target freq table...\n")
freqTable <- summarize_gene_pairs(gn_top200_list)

#cat("reading in TF top 200 targets ...\n")
#cellType_top200 = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl", paste0(cell_type, "_tf_top200_cor_intg_list.rds" )))
#cellType_top200_summary <- process_cellType_top200(cellType_top200, avg = TRUE)

# 4. gene-target integrated summary =================================================================================================== 
# =====================================================================================================================================
cat("Generating gene-target integrated summary...\n")

#tfs = names(TFs)
#intg_summary <- integrate_gene_pairs_summaries(gn_top200_list, tfs, freqTable, condition, names(ribo_list), matched_data_mm)
intg_summary <- integrate_gene_pairs_summaries(gn_allTrg_list, names(TFs), freqTable, condition, names(ribo_list), matched_data_mm)

output_file3 <- file.path(output_dir, paste0(cell_type, "_", condition, "_", mm_agg, "_hg_mm_intg_summary.rds"))
cat("Saving hg mm intg. summary to:", output_file3, "\n")
saveRDS(intg_summary, output_file3)

# 5. gene hg mm s.cor =================================================================================================== 
# =====================================================================================================================================
cat("Generating geneA hg mm s.cor summary...\n")
geneA_summary <- compute_geneA_summary(intg_summary)
geneA_summary <- geneA_summary %>% 
  mutate(is_TF = geneA %in% names(TFs)) %>%
  mutate(is_ribo = geneA %in% names(ribo_list))

output_file4 <- file.path(output_dir, paste0(cell_type, "_", condition, "_", mm_agg, "_geneA_hg_mm_sCor_summary.rds"))
cat("Saving geneA hg mm s.cor summary to:", output_file4, "\n")
saveRDS(geneA_summary, output_file4)


# Rscript 
# /cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list/all-matrix
# Mic
# CTL
# /home/nelazzabi/rewiring/data/coExpr/AD-Ctl 
# /home/nelazzabi/rewiring/data/coExpr/mouse/aggregate_cormat_allrank_microglia_mm.RDS
# mm_allrank
# /home/nelazzabi/rewiring/data/coExpr/mouse/aggregate_cormat_FZ_microglia_mm.RDS
# mm_FZ











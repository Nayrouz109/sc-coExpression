#!~/anaconda3/envs/biof-r/bin/R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(magrittr)
})

source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/extract_cor.r")



print("Reading Microglia data...")
mic_all_trgs <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_tf_allTrg_list.rds")
print("Reading Astrocyte data...")
ast_all_trgs <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Ast_tf_allTrg_list.rds")
print("Reading Excitatory Neuron data...")
exc_all_trgs <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Exc_tf_allTrg_list.rds")

names(mic_all_trgs) = gsub("^(Mic|Ast|Exc)_", "", names(mic_all_trgs))
names(ast_all_trgs) = gsub("^(Mic|Ast|Exc)_", "", names(ast_all_trgs))
names(exc_all_trgs) = gsub("^(Mic|Ast|Exc)_", "", names(exc_all_trgs))




result_df_test <- compute_signed_rank_norm_diff(mic_all_trgs, ast_all_trgs, exc_all_trgs)
saveRDS(result_df_test, "/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/mic_ast_exc_aggCorDiffV1.rds")

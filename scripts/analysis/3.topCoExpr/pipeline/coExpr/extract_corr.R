#!~/anaconda3/envs/biof-r/bin/R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(magrittr)
})


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript extract_corr.R <condition> <geneAA> <geneBB> <output_dir>")
}

condition <- args[1]
geneAA_input <- args[2]
geneBB_input <- args[3]
output_dir <- args[4]  

# Parse genes (assuming comma-separated input)
geneAA <- unlist(strsplit(geneAA_input, ","))
geneBB <- unlist(strsplit(geneBB_input, ","))



# source functions 
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/extract_cor.r")

# Load cell-type-specific lists
print("Reading Microglia data...")
#Mic <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_tf_allTrg_list.rds") #TF + all targets 
Mic <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_ribo_allGenes_allTrg_list.rds") # ribo + all targets 
Mic <- Mic[grep(paste0("_" , condition), names(Mic))]  
Mic <- bind_rows(Mic, .id = "dataset") 


print("Reading Astrocyte data...")
#Ast <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Ast_tf_allTrg_list.rds") #TF + all targets 
Ast <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Ast_ribo_allGenes_allTrg_list.rds") # ribo + all targets 
Ast <- Ast[grep(paste0("_" , condition), names(Ast))]  
Ast <- bind_rows(Ast, .id = "dataset") 


print("Reading Excitatory Neuron data...")
#Exc <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Exc_tf_allTrg_list.rds") #TF + all targets
Exc <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Exc_ribo_allGenes_allTrg_list.rds") # ribo + all targets 
Exc <- Exc[grep(paste0("_" , condition), names(Exc))] 
Exc <- bind_rows(Exc, .id = "dataset") 

print("Creating a list of cell types...")
cell_type_lists <- list(Mic = Mic, Ast = Ast, Exc = Exc)
rm(Mic)
rm(Ast)
rm(Exc)

#print("Reading catg data...")
#top_bottom_combined <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/mic_ast_exc_aggCorDiff_catg.rds") %>% {filter(.,.$category != "no_catg")}
print("Reading ribo data...")
ribo = fread("/home/nelazzabi/rewiring/data/genes/ribosomal/group-1054.csv", header = TRUE) %>% filter(Group %in% c("L ribosomal proteins", "S ribosomal proteins")) 
ribo_list <- split(ribo, ribo$`Approved symbol`)  
names(ribo_list) <- make.unique(names(ribo_list))


print("Extracting rank-normalized correlations...")
#final_result <- extract_rank_norm_corr(cell_type_lists, geneAA, geneBB)
#final_result <- extract_rank_norm_corrV1(cell_type_lists, top_bottom_combined)
final_result <- extract_rank_norm_corr_ribo(cell_type_lists, names(ribo_list))


# Save results to output directory
output_file <- file.path(output_dir, paste0("rank_corr_ribo", condition, ".csv"))
print(paste("Saving results to:", output_file))
fwrite(final_result, output_file)

print("Done!")


# Rscript 
# extract_corr.R
# CTL
# MEF2C
# PRKCA
# /home/nelazzabi/rewiring/data/coExpr/AD-Ctl/int_genes 

#Rscript extract_corr.R control "MEF2A,MEF2A,MEF2C,MEF2C,RORA,RORA,NFIB,RORB,RORB" "DOCK4,SFMBT2,PRKCA,ITPR2,ARHGAP24,SLC1A3,RORB,BCAR3,TGFBR1" /path/to/output/



library(tidyverse)
library(magrittr)
library(dplyr)
library(Matrix)
library(Seurat)


file <- '/cosmos/data/downloaded-data/AD-meta/ROSMAP-Micro/data/0-source/expr/ROSMAP.ImmuneCells.6regions.snRNAseq.counts.rds'
meta <- '/cosmos/data/downloaded-data/AD-meta/ROSMAP-Micro/data/0-source/expr/ROSMAP.ImmuneCells.6regions.snRNAseq.meta.rds'
out <- '/space/scratch/bxu_project/ROSMAP-Micro/'

print( paste("reading", basename(file))) 

data = readRDS(file) 
metadata <- readRDS(meta) 
mtxExtracted = data

metadata$patientID <- metadata$subject

metaExtracted <- as.data.frame(metadata) %>% 
  dplyr::arrange(subject) %>% mutate(patientID = factor(subject))

if (all( rownames(metaExtracted) %in% colnames(mtxExtracted))) {
  mtxExtracted <- mtxExtracted[, as.character(rownames(metaExtracted))]
} else {
  stop("Mismatch between patientIDs in metaExtracted and column names in mtxExtracted")
}


name   = paste0("ROSMAP-Micro") 

if (all(rownames(metaExtracted) == colnames(mtxExtracted))){
  print( paste("saving", basename(file))) 
  
  Matrix::writeMM(mtxExtracted, paste0(out, name, "_counts.mtx"))
  write.csv(rownames(mtxExtracted), paste0(out, name, "_genes.csv"), row.names = FALSE)
  write.csv(colnames(mtxExtracted), paste0(out, name, "_cell_ids.csv"), row.names = FALSE)
  write.csv(metaExtracted, paste0(out, name, "_meta_data.csv"), row.names = TRUE)
}
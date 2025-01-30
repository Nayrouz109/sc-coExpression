

#!~/anaconda3/bin/Rscript

library(tidyverse)
library(magrittr)
library(dplyr)
library(Matrix)
library(Seurat)


dir = "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/expr" 
out = "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/expr/rds-comps/"

files <- list.files(path = dir, pattern = "\\.rds$", full.names = TRUE) 

for (file in files) {
    print( paste("reading", basename(file))) 

    data = readRDS(file) 
    mtxExtracted = data$expr
    data$meta$patientID = data$meta$projid
    metaExtracted <- as.data.frame(data$meta) %>% 
                dplyr::arrange(patientID) %>% mutate(patientID = factor(patientID))
    
    if (all( rownames(metaExtracted) %in% colnames(mtxExtracted))) {
      mtxExtracted <- mtxExtracted[, as.character(rownames(metaExtracted))]
    } else {
      stop("Mismatch between patientIDs in metaExtracted and column names in mtxExtracted")
    }

    name0  = gsub("_.*" , "", basename(file)) 
    name   = paste0(name0, "_ROSMAP-Mathys-2019") 

    if (all(rownames(metaExtracted) == colnames(mtxExtracted))){
        print( paste("saving", basename(file))) 
        
        Matrix::writeMM(mtxExtracted, paste0(out, name, "_counts.mtx"))
        write.csv(rownames(mtxExtracted), paste0(out, name, "_genes.csv"), row.names = FALSE)
        write.csv(colnames(mtxExtracted), paste0(out, name, "_cell_ids.csv"), row.names = FALSE)
        write.csv(metaExtracted, paste0(out, name, "_meta_data.csv"), row.names = TRUE)
    }
}


# run like this 
# /space/opt/bin/Rscript extract-rds-kamath.r
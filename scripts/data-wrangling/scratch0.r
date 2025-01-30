
#!~/anaconda3/envs/testNE/bin/Rscript

### 00README 
#1. this script is to make sure files that have "MyelinRemoved" have the same genes as their cell-type corresponding files 
#2. integrates "MyelinRemoved" with the other ling files (per cell type)

#3. removes genes that have 0 in all cells (control and disease) per cell type 
#4. filters out genes that are detected in less than 0.02% of all single cells 
#5. filter out single cells that express >0.02% of genes  

# Load necessary libraries
library(Matrix)
library(data.table)
library(dplyr)
library(parallel)  # For mclapply



#### Functions ==============================================================================================
source("~/benchmark/co-expression/mathys/functions/filterSCmatrix-simple.R")

update_expr_matrix <- function(total_genes, expr_matrix) {
  
  missing_genes <- setdiff(total_genes, rownames(expr_matrix))
  if (length(missing_genes) > 0) {
    missing_rows <- matrix(0, nrow = length(missing_genes), ncol = ncol(expr_matrix))
    colnames(missing_rows) <- colnames(expr_matrix)
    rownames(missing_rows) <- missing_genes
    
    expr_matrix <- rbind(expr_matrix, missing_rows) # Bind the missing rows to the expression matrix
    
  }
  
  expr_matrix <- expr_matrix[rownames(expr_matrix) %in% total_genes, ]
  
  # Sanity check
  if (!(length(total_genes) == length(rownames(expr_matrix)) && all(total_genes %in% rownames(expr_matrix)))) {
    stop("Error: Length of total_genes doesn't match rownames of expr_matrix or some genes are missing.")
  }
  # Return the updated expression matrix
  return(expr_matrix)
}


###### Setup ==============================================================================
cell_types <- c("gabaergic", "glutamatergic", "astrocyte", "oligodendrocyte", "microglia", "polydendrocyte")
dir <- "/space/scratch/nairuz-rewiring-normal/data/processed/ling"


#MyelinFile_paths <- list.files(path = dir, pattern = "*MyelinRemoved.rds", full.names = TRUE)
allFile_paths <- list.files(path = dir, pattern = "*.rds", full.names = TRUE)






###### integrate myline with other ling ==============================================================================


for (cellType in cell_types) {
    
    files <- list.files(path = dir, pattern = paste0(cellType, ".*\\.rds$"), full.names = TRUE) 
    all_genes <- Reduce(union, lapply(files, function(file) rownames(readRDS(file)$expr)))

    if (cellType == "glutamatergic") {


    } else {
        updated_files_list <- lapply(files, function(file) {
            file_expr <- readRDS(file)$expr
            file_expr_updated <- update_expr_matrix(all_genes, file_expr)
            return(file_expr_updated)
        })
        
        combined_expr_matrix <- do.call(rbind, updated_files_list)
        combined_expr_matrixFilt <-  cleanCtmat(combined_expr_matrix, geneThr = 0.02, sampleThr = 0.02)
    } 
}



### fix meta data 
### add print statement 
### save the updated files somewhere new 

cellType = cell_types[2]


files <- list.files(path = dir, pattern = paste0(cellType, ".*\\.rds$"), full.names = TRUE) 
all_genes <- Reduce(union, lapply(files, function(file) rownames(readRDS(file)$expr)))

updated_files_list <- lapply(files, function(file) {
  file <- readRDS(file)
  file_expr_updated = list()
  file_expr_updated$expr <- update_expr_matrix(all_genes, file$expr)
  file_expr_updated$meta <- file$meta
  return(file_expr_updated)
})

glut1 <- updated_files_list[[1]]
glut2 <- updated_files_list[[2]]
glut3 <- updated_files_list[[3]]
glut4 <- updated_files_list[[4]]


mix1 = list(meta = rbind(glut1$meta, glut2$meta) , expr = cbind(glut1$expr, glut2$expr))
mix2 = list(meta = rbind(glut3$meta, glut4$meta) , expr = cbind(glut3$expr, glut4$expr))


mix1$expr <- cleanCtmat(mix1$expr, geneThr = 0.02, sampleThr = 0.02)
mix2$expr <- cleanCtmat(mix2$expr, geneThr = 0.02, sampleThr = 0.02)

all_genes2 <- union(rownames(mix1$expr), rownames(mix2$expr)) 
mix1$expr  <- update_expr_matrix(all_genes2, mix1$expr)
mix2$expr  <- update_expr_matrix(all_genes2, mix2$expr)

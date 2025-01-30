#!~/anaconda3/envs/testNE/bin/Rscript

### 00README 
#1. This script ensures that files with "MyelinRemoved" have the same genes as their corresponding cell-type files.
#2. Integrates "MyelinRemoved" with other ling files (per cell type).
#3. Removes genes that have 0 expression across all cells (control and disease) per cell type.

# Load necessary libraries
library(Matrix)
library(data.table)
library(dplyr)
library(parallel)  # For mclapply

#### Functions ==============================================================================================

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
output_dir <- "/space/scratch/nairuz-rewiring/data/0.processed-raw"

cell_type_map <- list(
  gabaergic = "Inh",
  glutamatergic = "Exc",
  astrocyte = "Ast",
  oligodendrocyte = "Oli",
  microglia = "Mic",
  polydendrocyte = "Opc"
)

###### Integrate Myelin with other ling files ==============================================================================

for (cellType in cell_types) {
    
    print(paste("Processing cell type:", cellType))
    
    # Read in files and get all unique genes
    files <- list.files(path = dir, pattern = paste0(cellType, ".*\\.rds$"), full.names = TRUE)
    all_genes <- Reduce(union, lapply(files, function(file) rownames(readRDS(file)$expr)))
    
    message(paste("updating genes, files for:", cellType))
    updated_files_list <- lapply(files, function(file) {
        file <- readRDS(file)
        file_expr_updated = list()
        file_expr_updated$expr <- update_expr_matrix(all_genes, file$expr)
        file_expr_updated$meta <- file$meta
        return(file_expr_updated)
    })

    cell_type_abbr <- cell_type_map[[cellType]]
    filename <- file.path(output_dir, paste0(cell_type_abbr, "_ling-2024.rds"))

    if (cellType == "glutamatergic") {
        message("Processing glutamatergic batches")
        
        glut1 <- updated_files_list[[1]]
        glut2 <- updated_files_list[[2]]
        glut3 <- updated_files_list[[3]]
        glut4 <- updated_files_list[[4]]
        
        # Combine and clean first batch
        mix1 <- list(
          meta = rbind(glut1$meta, glut2$meta) %>% tibble::column_to_rownames("barcodes"), 
          expr = cbind(glut1$expr, glut2$expr) %>% { .[rowSums(.) != 0, ] }
        )
        
        # Combine and clean second batch
        mix2 <- list(
          meta = rbind(glut3$meta, glut4$meta) %>% tibble::column_to_rownames("barcodes"), 
          expr = cbind(glut3$expr, glut4$expr) %>% { .[rowSums(.) != 0, ] }
        )
        
        # Get all genes across both batches
        allgenes <- union(rownames(mix1$expr), rownames(mix2$expr))
        
        mix1V2 <- list(meta = mix1$meta, expr = update_expr_matrix(allgenes, mix1$expr))
        mix2V2 <- list(meta = mix2$meta, expr = update_expr_matrix(allgenes, mix2$expr))

        mix1V2$meta$patientID <- mix1V2$meta$DONOR
        mix2V2$meta$patientID <- mix2V2$meta$DONOR
        
        # Check if the column names of expr match row names of meta
        if (!all(colnames(mix1V2$expr) == rownames(mix1V2$meta))) {
            stop("Error: Column names of mix1V2 expr do not match row names of meta.")
        } else {
            print("mix1V2: Column names of expr match row names of meta.")
        }
        
        if (!all(colnames(mix2V2$expr) == rownames(mix2V2$meta))) {
            stop("Error: Column names of mix2V2 expr do not match row names of meta.")
        } else {
            print("mix2V2: Column names of expr match row names of meta.")
        }
        
        # Save the results
        filename1 <- file.path(output_dir, paste0(cell_type_abbr, "_batch1", "_ling-2024.rds"))
        filename2 <- file.path(output_dir, paste0(cell_type_abbr, "_batch2", "_ling-2024.rds"))
        message(paste("saving files for:", cellType))
        saveRDS(mix1V2, filename1)
        saveRDS(mix2V2, filename2)

    } else {
        # Extract meta and expr separately
        meta_list <- lapply(updated_files_list, function(x) x$meta)
        expr_list <- lapply(updated_files_list, function(x) x$expr)

        combdinedData <- list(
          meta = do.call(rbind, meta_list) %>% tibble::column_to_rownames("barcodes"), 
          expr = do.call(cbind, expr_list) %>% { .[rowSums(.) != 0, ] }
        )
        combdinedData$meta$patientID <- combdinedData$meta$DONOR

        
        # Check if the column names of expr match row names of meta
        if (!all(colnames(combdinedData$expr) == rownames(combdinedData$meta))) {
            stop("Error: Column names of combined expr do not match row names of meta.")
        } else {
            print("combinedData: Column names of expr match row names of meta.")
        }
        
        # Save the results
        message(paste("saving files for:", cellType))
        saveRDS(combdinedData, filename)
    }
}

print("Processing complete!")

#!~/anaconda3/envs/testNE/bin/Rscript

# Load necessary libraries
library(rhdf5)
library(Matrix)
library(data.table)
library(dplyr)
library(parallel)  # For mclapply

#### Command-line argument parsing ####
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide a pattern as an argument, e.g., '*MyelinRemoved.h5'")
}
pattern <- args[1]

#### functions ==============================================================================================
process_h5_file <- function(h5file, cell_barcodes, cell_meta, prefix) {
  # Function to read an h5 file and extract barcodes and gene expression values

  barcodes <- h5read(h5file, "/matrix/barcodes")
  data <- h5read(h5file, "/matrix/data")
  indices <- h5read(h5file, "/matrix/indices")
  indptr <- h5read(h5file, "/matrix/indptr")
  shape <- h5read(h5file, "/matrix/shape")
  gene_names <- h5read(h5file, "/matrix/features/name")

  # Find the relevant barcodes that are present in the test1 dataframe
  relevant_cells <- barcodes[(barcodes %in% cell_barcodes)]
  
  if (length(relevant_cells) == 0) return(NULL)

  # Create a sparse matrix for the relevant cells
  expr_matrix <- sparseMatrix(i = indices + 1, p = indptr, x = data, dims = c(shape[1], shape[2]))
  colnames(expr_matrix) <- barcodes
  rownames(expr_matrix) <- gene_names

  # Subset to the relevant expression data (for the relevant cells)
  relevant_expr <- expr_matrix[, colnames(expr_matrix) %in% relevant_cells]
  colnames(relevant_expr) <- colnames(expr_matrix)[colnames(expr_matrix) %in% relevant_cells]
  
  # Create a meta dataframe for these cells based on test1_df
  meta <- cell_meta %>% filter(CELL_BARCODE %in% relevant_cells | barcodes %in% relevant_cells) %>% mutate(file_prefix = prefix)

  # Return list of meta and expr
  return(list(meta = meta, expr = relevant_expr))
}

align_to_common_genes <- function(expr_matrix, all_genes, file_genes) {
  # Function to align matrices to a common gene set
  aligned_matrix <- Matrix(0, nrow = length(all_genes), ncol = ncol(expr_matrix), sparse = TRUE)
  rownames(aligned_matrix) <- all_genes

  # Find the matching indices of the file's genes in the full gene list
  matching_indices <- match(file_genes, all_genes)
  
  # Insert the expression data into the correct rows
  aligned_matrix[matching_indices, ] <- expr_matrix
  colnames(aligned_matrix) <- colnames(expr_matrix)
  
  return(aligned_matrix)
}

###### getting ready ==============================================================================
# Define directory and file paths
#cell_types <- c("gabaergic", "glutamatergic", "astrocyte", "endothelia", "oligodendrocyte", "microglia", "polydendrocyte")
#cell_types <- c("astrocyte", "endothelia", "oligodendrocyte", "microglia", "polydendrocyte") # first one - done 
cell_types <- c("gabaergic")
#cell_types <- c("gabaergic", "glutamatergic",
dir <- "/space/scratch/nairuz-rewiring-normal/data/raw/ling/counts/counts-messy"

# Get all files matching the provided pattern
file_paths <- list.files(path = dir, pattern = pattern, full.names = TRUE)
if (length(file_paths) == 0) stop("No files found matching the pattern")

# Read the metadata
cells_meta_all <- fread("/space/scratch/nairuz-rewiring-normal/data/raw/ling/BA46.jointMetadata.txt")
cells_meta_all <- cells_meta_all %>% mutate(barcodes = paste0(PREFIX, "_", CELL_BARCODE))

# Initialize lists for storing matrices per cell type
cell_type_matrices <- lapply(cell_types, function(x) list(meta = NULL, expr = NULL))

####### wrangle ==========================================================================================

# Step 1: Extract gene lists from all .h5 files
all_genes <- c()  # To store the union of all genes across files
file_genes_list <- list()  # Store genes for each file

for (file in file_paths) {
  file_genes <- h5read(file, "/matrix/features/name")
  file_genes_list[[basename(file)]] <- file_genes
  all_genes <- union(all_genes, file_genes)
}

# Step 2: Process each .h5 file and align matrices using mclapply for parallelization
num_cores <- detectCores() - 1  # Use all but one core
process_file <- function(file) {
  prefix <- gsub(".h5", "", basename(file))
  message(paste("working on file", prefix))
  
  results <- lapply(cell_types, function(cell_type) {
    cell_type_meta <- cells_meta_all[predClass == cell_type & PREFIX %in% prefix]
    cell_type_barcodes <- cell_type_meta$CELL_BARCODE

    message(paste("working on cell type", cell_type))
    result <- process_h5_file(file, cell_type_barcodes, cell_type_meta, prefix)

    # Reorder meta rows and check barcode alignment
    result$meta <- result$meta[match(colnames(result$expr), result$meta$CELL_BARCODE), ]
    colnames(result$expr) <- result$meta$barcodes

    # Check if column names in expr match barcodes in meta
    if (!all(colnames(result$expr) == result$meta$barcodes)) {
      stop(paste("Error: Mismatch between expression matrix column names and meta barcodes in file", prefix, "for cell type", cell_type))
    }

    if (!is.null(result)) {
      file_genes <- file_genes_list[[basename(file)]]
      aligned_expr <- align_to_common_genes(result$expr, all_genes, file_genes)

      # Append meta and expression data to corresponding cell type
      list(meta = result$meta, expr = aligned_expr)
    } else {
      NULL
    }
  })
  
  return(results)
}

# Run the processing in parallel using mclapply
all_results <- mclapply(file_paths, process_file, mc.cores = num_cores)

# Step 3: Combine results per cell type
for (i in seq_along(cell_types)) {
  cell_type <- cell_types[[i]]
  
  # Combine meta and expr for each cell type from all files
  meta_combined <- do.call(rbind, lapply(all_results, function(res) res[[i]]$meta))
  expr_combined <- do.call(cbind, lapply(all_results, function(res) res[[i]]$expr))
  
  cell_type_matrices[[cell_type]]$meta <- meta_combined
  cell_type_matrices[[cell_type]]$expr <- expr_combined
}

# Step 4: Save matrices in RDS files
#output_dir <- "/space/scratch/nairuz-rewiring-normal/data/processed/ling"
output_dir <- "/space/scratch/nairuz-rewiring-normal/data/processed/ling/MyelinRemoved"
for (cell_type in cell_types) {
  output_file <- paste0(output_dir, "/", cell_type, "_ling.rds")
  saveRDS(cell_type_matrices[[cell_type]], output_file)
  message(paste("Saved", cell_type, "matrix to", output_file))
}



#### run the script 
#chmod +x ling-wrangling-parallel.r
#Rscript ling-wrangling-parallel.r "*MyelinRemoved.h5"
#Rscript ling-wrangling-parallel.r "BA46_(05-07-2019|05-21-2021|05-28-2021|07-10-2019|08-09-2019|08-27-2019|09-02-2021|09-03-2021|09-18-2019|09-22-2021|10-16-2019)_rxn[0-9]+\\.h5"


### sanity check to make sure total number of cell IDs in the big cell meta-data file matches the newly wrangled files 
# dim(cells_meta_all)
# 1228450      25
# length(unique(cells_meta_all$PREFIX))
# 87
# length(unique(cells_meta_all$CELL_BARCODE))
# 1048695
# length(unique(cells_meta_all$barcodes))
# 1,228,450

## newly wrangled with pattern BA46_(date)_rxn[0-9]+\\.h5"
# Cell Type: Astrocyte 
# Expr dimensions: 39735 177279 

# Cell Type: Endothelia 
# Expr dimensions: 39735 17762 

# Cell Type: GABAergic 
# Expr dimensions: 39735 224485 

# Cell Type: Microglia 
# Expr dimensions: 39735 43322 

# Cell Type: Oligodendrocyte 
# Expr dimensions: 39735 130205 

# Cell Type: Polydendrocyte 
# Expr dimensions: 39735 63607 
# 177279 + 63607 + 130205 + 43322 + 224485 + 17762 = 656,660 


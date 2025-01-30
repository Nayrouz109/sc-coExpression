

#!~/anaconda3/envs/testNE/bin/Rscript

# Load necessary libraries
library(rhdf5)
library(Matrix)
library(data.table)
library(dplyr)

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
  relevant_expr <- expr_matrix[, colnames(expr_matrix) %in% relevant_cells]
  
  # Create a meta dataframe for these cells based on test1_df
  meta <- cell_meta %>% filter(CELL_BARCODE %in% relevant_cells | barcodes %in% relevant_cells) %>% mutate(file_prefix = prefix)

  # Return list of meta and expr
  return(list(meta = meta, expr = relevant_expr))
}

align_to_common_genes <- function(expr_matrix, all_genes, file_genes) {
  # Function to align matrices to a common gene set
  aligned_matrix <- Matrix(0, nrow = length(all_genes), ncol = ncol(expr_matrix), sparse = TRUE)
  
  # Find the matching indices of the file's genes in the full gene list
  matching_indices <- match(file_genes, all_genes)
  
  # Insert the expression data into the correct rows
  aligned_matrix[matching_indices, ] <- expr_matrix
  
  return(aligned_matrix)
}

###### getting ready ==============================================================================
# Define directory and file paths
cell_types <- c("gabaergic", "glutamatergic", "astrocyte", "endothelia", "oligodendrocyte", "microglia", "polydendrocyte")
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

# Step 2: Process each .h5 file and align matrices
for (file in file_paths) {
  #prefix <- gsub("_MyelinRemoved.h5", "", basename(file))
  prefix <- gsub(".h5", "", basename(file))
  message(paste("working on file", prefix))
  
  for (cell_type in cell_types) {
    cell_type_meta <- cells_meta_all[predClass == cell_type & PREFIX %in% prefix]
    cell_type_barcodes <- cell_type_meta$CELL_BARCODE

    message(paste("working on cell type", cell_type))
    result <- process_h5_file(file, cell_type_barcodes, cell_type_meta, prefix)
    
    if (!is.null(result)) {
      file_genes <- file_genes_list[[basename(file)]]
      aligned_expr <- align_to_common_genes(result$expr, all_genes, file_genes)

      # Append meta and expression data to corresponding cell type
      cell_type_matrices[[cell_type]]$meta <- rbind(cell_type_matrices[[cell_type]]$meta, result$meta)
      cell_type_matrices[[cell_type]]$expr <- cbind(cell_type_matrices[[cell_type]]$expr, aligned_expr)
    }
  }
}

# Step 3: Save matrices in RDS files
output_dir <- "/space/scratch/nairuz-rewiring-normal/data/processed/ling"
for (cell_type in cell_types) {
  output_file <- paste0(output_dir, "/", cell_type, "_ling.rds")
  saveRDS(cell_type_matrices[[cell_type]], output_file)
  message(paste("Saved", cell_type, "matrix to", output_file))
}



#### run the script 
#chmod +x ling-wrangling.r
#Rscript ling-wrangling.r "*MyelinRemoved.h5"
#Rscript ling-wrangling.r "BA46_(05-07-2019|05-21-2021|05-28-2021|07-10-2019|08-09-2019|08-27-2019|09-02-2021|09-03-2021|09-18-2019|09-22-2021|10-16-2019)_rxn[0-9]+\\.h5"


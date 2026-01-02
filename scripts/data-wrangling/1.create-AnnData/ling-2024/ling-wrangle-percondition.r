#!~/anaconda3/envs/testNE/bin/Rscript

# Load necessary libraries
library(Matrix)
library(data.table)
library(dplyr)
library(parallel)  # For mclapply


#### functions ==============================================================================================
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

get_patient_category <- function(ptn_id, ptns_meta) {
    # categorize patients as control or schizophrenia based on metadata
  patient_row <- ptns_meta %>% filter(nbb.id == ptn_id)
  if (nrow(patient_row) == 0) {
    return(NA)
  }
  if (patient_row$Schizophrenia == "Affected") {
    return("schizophrenia")
  } else {
    return("control")
  }
}

###### Setup ==============================================================================
cell_types <- c("gabaergic", "glutamatergic", "astrocyte", "oligodendrocyte", "microglia", "polydendrocyte")
dir <- "/space/scratch/nairuz-rewiring-normal/data/processed/ling"
output_dir <- "/space/scratch/nairuz-rewiring/data/processed-raw"
file_paths <- list.files(path = dir, pattern = "*.rds" , full.names = TRUE)



all_genes <- c()  # To store the union of all genes across files
file_genes_list <- list()  # Store genes for each file


for (file in file_paths) {
    print(paste0("storing genes for:", basename(file))) 

    file_genes <- rownames(readRDS(file)$expr)  
    file_genes_list[[basename(file)]] <- file_genes
    all_genes <- union(all_genes, file_genes)
}

# Read the metadata
cells_meta_all <- fread("/space/scratch/nairuz-rewiring-normal/data/raw/ling/BA46.jointMetadata.txt")
cells_meta_all <- cells_meta_all %>% mutate(barcodes = paste0(PREFIX, "_", CELL_BARCODE))

ptns_meta_all <- fread("/space/scratch/nairuz-rewiring-normal/data/raw/ling/SZvillage_donorMetadata.txt")

cell_type_map <- list(
  gabaergic = "Inh",
  glutamatergic = "Exc",
  astrocyte = "Ast",
  oligodendrocyte = "Oli",
  microglia = "Mic",
  polydendrocyte = "Opc"
)





# Function to process files per cell type
process_celltype_files <- function(cell_type, meta_data, dir) {
  
  message(paste0("Processing files for cell type: ", cell_type))

  files <- list.files(path = dir, pattern = paste0(cell_type, "_ling"), full.names = TRUE)   # Find all files related to the current cell type
  
  if (length(files) == 0) {
    message(paste0("No files found for cell type: ", cell_type))
    return(NULL)
  } else {
    message(paste0(length(files), " files found for cell type: ", cell_type))
  }

  # Initialize lists for control and schizophrenia data
  control_list <- list(meta = NULL, expr = NULL)
  schizophrenia_list <- list(meta = NULL, expr = NULL)
  
  for (file in files) {

    message(paste0("Loading data from file: ", file))
    data <- readRDS(file) 
    file_expr <- data$expr

    file_genes <- file_genes_list[[basename(file)]]
    aligned_expr <- align_to_common_genes(file_expr, all_genes, file_genes)
    
    # Extract patient ID from metadata
    patient_ids <- unique(data$meta$DONOR)
    
    message("Classifying patients into control or schizophrenia...")
    # Classify patients and separate into control and schizophrenia
    for (i in 1:length(patient_ids)) {
        patient_id <- patient_ids[i]
        category <- get_patient_category(patient_id, meta_data)
        
        donor_meta <- data$meta[ data$meta$DONOR == patient_id , ] %>% tibble::column_to_rownames("barcodes")
        donor_meta$patientID <- donor_meta$DONOR
        donor_expr <- file_expr[, colnames(file_expr) %in% rownames(donor_meta)]


        if (!all(rownames(donor_meta) == colnames(donor_expr))) {
            stop(paste0("Error: Rownames of donor_meta do not match colnames of donor_expr for patient ID: ", patient_id))
        }
      
        if (category == "control") {
            control_list$meta <- rbind(control_list$meta, donor_meta)
            control_list$expr <- cbind(control_list$expr, donor_expr)

        } else if (category == "schizophrenia") {
            schizophrenia_list$meta <- rbind(schizophrenia_list$meta, donor_meta)
            schizophrenia_list$expr <- cbind(schizophrenia_list$expr, donor_expr)
        }
    }
  }
  
  # Save control and schizophrenia lists to .rds files
  cell_type_abbr <- cell_type_map[[cell_type]]
  control_filename <- file.path(output_dir, paste0(cell_type_abbr, "_ling-2024_control.rds"))
  schizophrenia_filename <- file.path(output_dir, paste0(cell_type_abbr, "_ling-2024_schizophrenia.rds"))
  
  message(paste0("Saving control data to: ", control_filename))
  saveRDS(control_list, control_filename)
  
  message(paste0("Saving schizophrenia data to: ", schizophrenia_filename))
  saveRDS(schizophrenia_list, schizophrenia_filename)

}


# Use mclapply to parallelize the processing of each cell type
num_cores <- 6  # one core per cell type 
message(paste0("Using ", num_cores, " cores for parallel processing."))

mclapply(cell_types, function(cell_type) {
  process_celltype_files(cell_type, ptns_meta_all, dir)
}, mc.cores = num_cores)



message("All cell types have been processed successfully.")

#chmod +x ling-wrangle-percondition.r
#Rscript ling-wrangle-percondition.r 








##### read in what was provided by the authors 
ROSMAP_DeJager2023cells <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv") 
ROSMAP_DeJager2023cells$patientID = ROSMAP_DeJager2023cells$specimenID
ROSMAP_DeJager2023genes <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/ROSMAP_Brain.snRNAseq_metadata_genes_20230420.csv")


#### wranlge the matrix format provided 
library(Matrix)
data <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/ROSMAP_Brain.snRNAseq_counts_sparse_format_20230420.csv.gz", fill=FALSE)
data <- fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/ROSMAP_Brain.snRNAseq_counts_sparse_format_20230420.csv.gz", 
              skip = 1, # Skip the first line
              header = FALSE, # No header since it's a sparse format
              fill = TRUE)

colnames(data) <- c("row", "col", "value")
# Convert columns to numeric
data$row <- as.numeric(data$row)
data$col <- as.numeric(data$col)
data$value <- as.numeric(data$value)

data <- data[complete.cases(data), ] # Remove any rows with NA (if they exist)

n_rows <- max(data$row, na.rm = TRUE) # Determine matrix dimensions (based on the maximum row and column indices)
n_cols <- max(data$col, na.rm = TRUE)

sparse_mat <- sparseMatrix(i = data$row, 
                           j = data$col, 
                           x = data$value, 
                           dims = c(n_rows, n_cols)) # Create a sparse matrix (dgCMatrix)


rownames(sparse_mat) <- ROSMAP_DeJager2023genes$V2[-1]
colnames(sparse_mat) <- ROSMAP_DeJager2023cells$cell_name

all(colnames(sparse_mat) == ROSMAP_DeJager2023cells$cell_name)
#### split the data into multiple matrices based on cell type annotations 
library(Matrix)
library(data.table)


cell_types <- unique(ROSMAP_DeJager2023cells$`Cell Type`)
cell_type_lists <- list()

# Loop over each cell type and create a list containing expression and meta data
for (cell_type in cell_types) {

  meta_subset <- ROSMAP_DeJager2023cells[`Cell Type` == cell_type]
  expr_subset <- sparse_mat[, colnames(sparse_mat) %in% meta_subset$cell_name]
  
  cell_type_list <- list(
    expr = expr_subset,
    meta = meta_subset
  )
  
  # Store the list in the main list using the cell type as the key
  cell_type_lists[[cell_type]] <- cell_type_list
}



output_dir <- "/space/scratch/nairuz-rewiring/data/0.processed-raw"
cell_type_mapping <- list(
  "Exc" = "Exc",
  "Inh" = "Inh",
  "opcs" = "Opc",
  "astrocytes" = "Ast",
  "oligodendrocytes" = "Oli",
  "microglia" = "Mic"
)


for (cell_type in names(cell_type_mapping)) {

  short_label <- cell_type_mapping[[cell_type]]
  filename <- paste0(output_dir, "/", short_label, "_ROSMAP-DeJager-2023.rds")
  
  print(filename)
  
  # Save the RDS file
  saveRDS(cell_type_lists[[cell_type]], filename)
}



### wranlge patinet meta data 
all_patients <- readRDS("/space/scratch/nairuz-rewiring/data/1.processed-5per/all-patient-metaData.rds")


ROSMAP_DeJager2023cells <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager-2023/data/0-source/expr/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv") 
ROSMAP_DeJager2023cells$patientID = ROSMAP_DeJager2023cells$specimenID

ROSMAP_DeJager2023_patients <- ROSMAP_DeJager2023cells [, c("patientID" , "specimenID")] %>% unique()

ROSMAP_DeJager2023_patients <- ROSMAP_DeJager2023_patients %>%
  mutate(
    sex = NA, 
    age = NA, 
    study = "ROSMAP-DeJager-2023", 
    brainRegion = "DLPFC",  
    #ADdiag2types = ifelse(  (grepl("Cdx1|Cog1", patientID) & grepl("pAD0|Path0", patientID)), "no", "yes"), 
    ADdiag2types = fifelse(
      grepl("Cdx1|Cog1", patientID) & grepl("pAD0|Path0", patientID), "no",
      fifelse(
        grepl("Cdx1|Cog1", patientID) & grepl("pAD1|Path1", patientID), "yes",
        fifelse(
          grepl("Cdx4|Cog4", patientID) & grepl("pAD1|Path1", patientID), "yes",
          fifelse(
            grepl("Cdx4|Cog4", patientID) & grepl("pAD0|Path0", patientID), "no", 
            NA_character_
          )
        )
      )
    ), 
    ADdiag3types = fifelse(
      grepl("Cdx1|Cog1", patientID) & grepl("pAD0|Path0", patientID), "control",
      fifelse(
        grepl("Cdx1|Cog1", patientID) & grepl("pAD1|Path1", patientID), "early",
        fifelse(
          grepl("Cdx4|Cog4", patientID) & grepl("pAD1|Path1", patientID), "late",
          fifelse(
            grepl("Cdx4|Cog4", patientID) & grepl("pAD0|Path0", patientID), "control", 
            NA_character_
          )
        )
      )
    )
  ) %>%
  select(c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study"))

all_patients <- all_patients %>% filter(study != "ROSMAP-DeJager-2023")

all_patients_updated <- rbind(all_patients, ROSMAP_DeJager2023_patients)

saveRDS(all_patients_updated, "/space/scratch/nairuz-rewiring/data/1.processed-5per/all-patient-metaData.rds")








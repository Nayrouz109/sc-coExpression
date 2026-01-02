
library(dplyr)
library(stringr)


test = read_rds("/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/eric_queries/rds/ramos_raw.rds")
cln_meta = read_csv("/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/ramos/clean_metadata.csv")

test@meta.data <- test@meta.data %>%
  mutate(
    sample_suffix = str_extract(sample, "[0-9]+[A-Z]+$"), # get "9C" from "GSM6720852_9C"
    cellids = paste0(orig.ident, "_", sample_suffix)       # concatenate
  )

test@meta.data$celltypes <- cln_meta$celltypes[
  match(test@meta.data$cellids, cln_meta$cellids)
]


test@meta.data <- test@meta.data[
  test@meta.data$cellids %in% cln_meta$cellids,]




library(Seurat)
# Extract counts matrix
counts_matrix <- Seurat::GetAssayData(test, assay = "RNA", slot = "counts")
# Subset columns (cells)
cells_to_keep <- rownames(test@meta.data)
counts_subset <- counts_matrix[, colnames(counts_matrix) %in% cells_to_keep]
# Create a new Assay and assign
test[["RNA"]] <- Seurat::CreateAssayObject(counts = counts_subset)












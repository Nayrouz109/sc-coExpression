





dir = "/space/scratch/nairuz-rewiring-normal/data/processed/ling/processedV2"
output_dir = "/cosmos/data/downloaded-data/AD-meta/ling-2024/data/extracted-rds-comps/"

files = list.files(path = dir, pattern = "\\.rds$", full.names = TRUE)

### Extract Components for annData ===================================================
for (file in files) {
  cell_type <- gsub(".*/|\\.rds$", "", file)
  data = readRDS(file)
  metaExtracted <- data$meta %>% mutate (patientID = DONOR) 
  mtxExtracted  <- data$expr 
  
  all(colnames(mtxExtracted) %in% rownames(metaExtracted))
  
  file_name = paste0(output_dir, cell_type)
  
  Matrix::writeMM(mtxExtracted, paste0(output_dir, cell_type, "_countsV2.mtx"))
  write.csv(rownames(mtxExtracted), paste0(output_dir, cell_type, "_genes.csv"), row.names = FALSE)
  write.csv(colnames(mtxExtracted), paste0(output_dir, cell_type, "_cell_ids.csv"), row.names = FALSE)
  write.csv(metaExtracted, paste0(output_dir, cell_type, "_meta_data.csv"), row.names = TRUE)
}



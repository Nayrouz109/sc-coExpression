
#!~/anaconda3/envs/testNE/bin/R

#library(tidyverse)
#library(magrittr)
#library(dplyr)

#file_directory <- "/space/scratch/nairuz-rewiring/data/0.processed-raw"
file_directory <- "/space/scratch/nairuz-rewiring-normal/data/processed/ling/processedV2"
rds_files <- list.files(path = file_directory, pattern = "\\.rds$", full.names = TRUE)

missing_patientID_files <- list()


for (file in rds_files) {
  print(paste( "reading" , file ))

  data <- readRDS(file)
  

  if (!("patientID" %in% colnames(data$expr))) {
    missing_patientID_files <- append(missing_patientID_files, file)
  }
}

# Print the list of files missing "patientID" in data$expr
missing_patientID_files

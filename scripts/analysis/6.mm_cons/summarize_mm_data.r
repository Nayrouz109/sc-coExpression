#!~/anaconda3/envs/biof-r/bin/R

library(dplyr)
input_dir <- "/cosmos/data/project-data/NW-rewiring/data/mouse/0.original"

# Get list of subdirectories (studies)
study_dirs <- list.dirs(input_dir, recursive = FALSE)
num_studies <- length(study_dirs)  # Total number of studies

# Initialize an empty data frame to store the summary
summary_table <- data.frame(Study = character(), Num_Cells = integer(), stringsAsFactors = FALSE)

# Loop over each study directory
# Loop over each study directory
for (i in seq_along(study_dirs)) {
  study_dir <- study_dirs[i]
  study_name <- basename(study_dir)  # Extract study name
  
  # Find the target RDS file
  rds_file <- list.files(study_dir, pattern = "_clean_mat_and_meta.RDS$", full.names = TRUE)
  
  # Print processing status
  cat(sprintf("Processing %d/%d: %s\n", i, num_studies, study_name))
  
  # Ensure there is exactly one matching file before proceeding
  if (length(rds_file) == 1) {
    data <- readRDS(rds_file)
    
    # Check if Mat exists and has correct dimensions
    if ("Mat" %in% names(data)) {
      num_cells <- dim(data$Mat)[2]  # Extract number of single cells
      print(num_cells)
      
      # Append to summary table
      summary_table <- rbind(summary_table, data.frame(Study = study_name, Num_Cells = num_cells))
    } else {
      warning(sprintf("Skipping %s: 'Mat' is missing or not a matrix.\n", study_name))
    }
  } else {
    warning(sprintf("Skipping %s: No valid RDS file found.\n", study_name))
  }
}



# Save the summary table using the input directory
output_file <- file.path(input_dir, "mm_study_summary.csv")
write.csv(summary_table, file = output_file, row.names = FALSE)



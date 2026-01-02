
#!~/anaconda3/envs/biof-r/bin/R

# Load required packages
suppressMessages({
  library(Matrix)
  library(tidyverse)
  library(stringr)
})


# Define the input directory
data_dir <- "/cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/PB/cpm-log"
cat("📁 Using data directory:", data_dir, "\n")

# List .rds files that contain "Rex"
rds_files <- list.files(data_dir, full.names = TRUE)
cat("📄 Found", length(rds_files), "RDS files.\n")

# Example tf_df definition — REPLACE with actual file or source
# You must replace this with how you load or define tf_df in your project
tf_df = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/TF_rpr_inegrated_cellTypes.rds")

# Initialize list to collect results
results_list <- list()

# Loop through expression files
for (file in rds_files) {
  cat("📥 Processing:", basename(file), "...\n")
  
  obj <- readRDS(file)
  expr <- obj$expr
  meta <- obj$meta
  
  # Parse metadata from filename
  file_base <- basename(file)
  cell_type <- str_extract(file_base, "^[^_]+")
  study_name <- str_extract(file_base, "(?<=_)[^\\.]+")
  
  # Subset genes
  genes_to_check <- intersect(tf_df$geneA, rownames(expr))
  
  if (length(genes_to_check) == 0) {
    cat("⚠️  No matching genes found in", file_base, "\n")
    next
  }
  expr = expr[ , colnames(expr) %in% rownames(meta %>% filter(disease == "CTL")) ]
  avg_expr <- Matrix::rowMeans(expr[genes_to_check, , drop = FALSE])
  
  df <- tibble(
    geneA = names(avg_expr),
    avg_expr_log_cpm = avg_expr,
    cell_type = cell_type,
    study_name = study_name
  )
  
  results_list[[file_base]] <- df
  cat("✅ Finished:", file_base, "\n")
}

# Combine results
final_expr_df <- bind_rows(results_list)
cat("📊 Final table contains", nrow(final_expr_df), "rows.\n")

# Optional: write to CSV or RDS
output_csv <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/TF_exp_summary_logCPM.csv"
write_csv(final_expr_df, output_csv)
cat("💾 Results written to", output_csv, "\n")

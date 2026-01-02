#!/home/nelazzabi/anaconda3/envs/biof-r/bin/R

# Load required packages
suppressMessages({
  library(Matrix)
  library(tidyverse)
  library(stringr)
})

# === Input Directory ===
data_dir <- "/cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/PB/cpm-log"
cat("📁 Using data directory:", data_dir, "\n")

# List RDS expression files
rds_files <- list.files(data_dir, full.names = TRUE, pattern = "\\.rds$")
cat("📄 Found", length(rds_files), "RDS files.\n")

# Initialize results list
results_list <- list()

# === Main Loop ===
for (file in rds_files) {
  cat("📥 Processing:", basename(file), "...\n")
  
  obj <- readRDS(file)
  expr <- obj$expr  # dgCMatrix
  meta <- obj$meta  # cell-level metadata

  # Parse metadata from filename
  file_base <- basename(file)
  cell_type <- str_extract(file_base, "^[^_]+")
  study_name <- str_extract(file_base, "(?<=_)[^\\.]+")

  # Subset to CTL cells
  ctl_cells <- rownames(meta %>% filter(disease == "CTL"))
  expr_ctl <- expr[, colnames(expr) %in% ctl_cells]

  # Skip if no CTL cells
  if (ncol(expr_ctl) == 0) {
    cat("⚠️  No CTL cells in", file_base, "\n")
    next
  }

  # Compute average expression for all genes
  avg_expr <- Matrix::rowMeans(expr_ctl)

  # Store results
  df <- tibble(
    gene = names(avg_expr),
    avg_expr_log_cpm = avg_expr,
    cell_type = cell_type,
    study_name = study_name, 
    condition = "CTL"
  )

  results_list[[file_base]] <- df
  cat("✅ Finished:", file_base, "with", nrow(df), "genes.\n")
}

# Combine into final table
final_expr_df <- bind_rows(results_list)
cat("📊 Final table contains", nrow(final_expr_df), "rows and", 
    length(unique(final_expr_df$gene)), "unique genes.\n")

# === Save Output ===
output_csv <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/allGenes_exp_summary_logCPM.csv"
#              "/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/allGenes_exp_summary_logCPM.csv"
write_csv(final_expr_df, output_csv)
cat("💾 Results written to", output_csv, "\n")

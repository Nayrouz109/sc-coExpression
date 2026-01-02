#!/home/nelazzabi/anaconda3/envs/biof-r/bin/R

# ────────────────────────────────────────────────────────────────────────────────
# TF Co-Expression & Unibind Integration Script
# ────────────────────────────────────────────────────────────────────────────────

suppressMessages({
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(tibble)
})


pwd <- "/home/nelazzabi/rewiring/data/coExpr/0all-human"
top200_files <- list.files(pwd, full.names = TRUE, pattern = "_tf_top200_list\\.rds$")
cat("📁 Using", length(top200_files), "input co-expression files.\n")

unibind_rds <- "/cosmos/data/project-data/NW-rewiring/data/unibind/09_Scored_Unibind_ChIPseq_data.RDS"

output_dir <- "/home/nelazzabi/rewiring/data/coExpr/0all-human-unibind"




# === Functions ===
get_cell_type <- function(path) str_extract(basename(path), "^[A-Za-z]+")

# ────────────────────────────────────────────────────────────────────────────────
# STEP 1: Summarize CTL pair frequencies per cell type
# ────────────────────────────────────────────────────────────────────────────────

CTL_freq_list <- map(top200_files, function(f) {
  cell_type <- get_cell_type(f)
  datalist <- readRDS(f)
  
  ctl_freq <- datalist %>%
    map_dfr(~ .x %>%
              filter(condition == "control") %>%
              select(geneA, geneB)
    ) %>%
    group_by(geneA, geneB) %>%
    summarise(freq = n(), .groups = "drop") %>%
    mutate(cell_type = cell_type)
  
  return(ctl_freq)
}) %>%
  set_names(map_chr(top200_files, get_cell_type))

# Combine all cell types
CTL_freq_all <- bind_rows(CTL_freq_list)

# ────────────────────────────────────────────────────────────────────────────────
# STEP 2: Pivot to wide format and compute pairwise differences
# ────────────────────────────────────────────────────────────────────────────────

CTL_freq_wide <- CTL_freq_all %>%
  unite("gene_pair", geneA, geneB, sep = "_", remove = FALSE) %>%
  select(gene_pair, geneA, geneB, cell_type, freq) %>%
  pivot_wider(names_from = cell_type, values_from = freq, values_fill = 0)

# Optional: compute pairwise differences between cell types
cell_types <- names(CTL_freq_list)
combs <- combn(cell_types, 2, simplify = FALSE)

for (combo in combs) {
  col1 <- combo[1]
  col2 <- combo[2]
  diff_col <- paste0("diff_", col1, "_", col2)
  CTL_freq_wide[[diff_col]] <- CTL_freq_wide[[col1]] - CTL_freq_wide[[col2]]
}

# Keep core summary for integration
TF_top_allCells <- CTL_freq_wide[, 1:9]

# ────────────────────────────────────────────────────────────────────────────────
# STEP 3: Process Unibind binding matrices
# ────────────────────────────────────────────────────────────────────────────────

cat("📦 Loading Unibind data from:", unibind_rds, "\n")
test <- readRDS(unibind_rds)
mat <- test$Permissive_hg$Mat_qnl
meta <- test$Permissive_hg$Meta
stopifnot(all(meta$ID %in% colnames(mat)))

# Map TF symbols to experiment columns
tf_by_col <- meta$Symbol[match(colnames(mat), meta$ID)]
tf_levels <- unique(tf_by_col)

# Compute average binding per TF
mat_by_tf <- sapply(tf_levels, function(tf) {
  cols <- which(tf_by_col == tf)
  if (length(cols) > 1) {
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  } else {
    mat[, cols]
  }
})

# Convert to data frame
mat_by_tf_df <- as.data.frame(mat_by_tf) %>% rownames_to_column("geneB")

# Create ranked version
mat_by_tf_ranked <- mat_by_tf_df %>%
  mutate(across(-geneB, ~ rank(.x, ties.method = "min")))

mat_by_tf_ranked_scaled <- mat_by_tf_df %>%
  mutate(across(
    -geneB,
    ~ {
      r <- rank(.x, ties.method = "min")
      (r - min(r)) / (max(r) - min(r))  # scale to [0,1]
    }
  ))

# ────────────────────────────────────────────────────────────────────────────────
# STEP 4: Integrative summaries (raw + ranked)
# ────────────────────────────────────────────────────────────────────────────────

# Function to join Unibind matrix (raw or ranked) with coexpression summary
integrate_with_binding <- function(mat_df, label) {
  mat_long <- mat_df %>%
    pivot_longer(cols = -geneB, names_to = "geneA", values_to = "binding_score") %>%
    mutate(gene_pair = paste(geneA, geneB, sep = "_"))
  
  merged <- TF_top_allCells %>%
    inner_join(mat_long, by = "gene_pair") %>%
    select(gene_pair, geneA = geneA.x, geneB = geneB.x, Ast, Exc, Inh, Mic, Oli, Opc, binding_score)
  
  cat("✅ Created integrative summary with", label, "binding scores: ",
      nrow(merged), "rows.\n")
  
  return(merged)
}

TF_top_with_binding_raw <- integrate_with_binding(mat_by_tf_df, "raw")
TF_top_with_binding_ranked <- integrate_with_binding(mat_by_tf_ranked, "ranked")
TF_top_with_binding_ranked_scaled <- integrate_with_binding(mat_by_tf_ranked_scaled, "ranked_scaled")
# ────────────────────────────────────────────────────────────────────────────────
# STEP 5: Save Output
# ────────────────────────────────────────────────────────────────────────────────

output_path <- file.path(output_dir, "TF_top_with_binding_ranked_scaled.rds")
saveRDS(TF_top_with_binding_ranked_scaled, output_path)
cat("💾 Saved ranked integrative summary to:", output_path, "\n")

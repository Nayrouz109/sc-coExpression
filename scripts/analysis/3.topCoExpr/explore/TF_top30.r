


dir <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell"
files <- list.files(dir, pattern = ".*rpr.*.rds", full.names = TRUE)
data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))



#####   extracting the top 30 TFs =============================================================
tf_elements <- data_list[grepl("t200_rpr_tf", names(data_list))]
top30_per_celltype <- lapply(tf_elements, function(df) {
  df %>%
    arrange(desc(meanTop200_AD)) %>%
    slice_head(n = 30)
})


top30_with_celltype <- Map(function(df, name) {
  df$cell_type <- sub("_t200_rpr_tf", "", name)
  df
}, top30_per_celltype, names(top30_per_celltype))


top30_summary_df <- bind_rows(top30_with_celltype)
write.csv(top30_summary_df, "/home/nelazzabi/temp/top30TF.csv")

#### which ones to ask Unibind to provide data for ==============================================================
tf_by_gene_binding <- readRDS("/home/nelazzabi/rewiring/data/modalities/gene_by_tf_binding_ranked_scaled.rds")
needed = (unique(top30_summary_df$geneA))[!(unique(top30_summary_df$geneA) %in% colnames(tf_by_gene_binding))]
needed_df = (top30_summary_df %>% filter(geneA %in% needed))

needed_summary_genes <- needed_df %>%
  group_by(geneA) %>%
  summarise(
    mic_TF   = any(mic_TF),              # TRUE if gene is TF in any row
    brain_TF = any(brain_TF),            # TRUE if brain TF in any row
    cell_type = paste(unique(cell_type), collapse = ", "),  # list all unique cell types
    top_freq  = n()                      # how many times gene appears
  ) %>%
  select(-c("mic_TF")) %>% 
  ungroup()

write.csv(needed_summary_genes, "/home/nelazzabi/rewiring/data/TFs/needed/topTF_missingUnibindData.csv")






##### cell type occurance ============================================================================================================================
# Add a column indicating if the gene is a TF
top30_summary_df <- top30_summary_df %>%
  mutate(is_TF = mic_TF | brain_TF)

# Filter only TFs
tf_df <- top30_summary_df %>% filter(is_TF)

# Count in how many cell types each TF appears
tf_counts <- tf_df %>%
  group_by(geneA) %>%
  summarise(cell_type_count = n_distinct(cell_type), .groups = "drop")

# Keep only TFs that appear in a single cell type
unique_tfs <- tf_counts %>% filter(cell_type_count == 1)

# Get the cell type for each unique TF and summarize per cell type
unique_tfs_summary <- tf_df %>%
  filter(geneA %in% unique_tfs$geneA) %>%
  group_by(cell_type) %>%
  summarise(unique_TFs = paste(geneA, collapse = ", "), .groups = "drop")

unique_tfs_summary




### 2 in inh and exc 
# Add a column indicating if the gene is a TF
top30_summary_df <- top30_summary_df %>%
  mutate(is_TF = mic_TF | brain_TF)

# Filter only TFs
tf_df <- top30_summary_df %>% filter(is_TF)

# Count in which cell types each TF appears
tf_cell_counts <- tf_df %>%
  #filter(cell_type %in% c("Inh", "Exc")) %>%
  group_by(geneA) %>%
  summarise(cell_types = paste(sort(unique(cell_type)), collapse = ", "),
            n_types = n_distinct(cell_type), .groups = "drop")

# Keep only TFs that appear in both Inh and Exc
tfs_in_both <- tf_cell_counts %>% filter(n_types == 2)

tfs_in_both



#####   average null per cell type  =============================================================

tn_names <- grep("t200_rpr_tn$", names(data_list), value = TRUE)
# Initialize empty list to store results
avg_overlap_summary <- data.frame(
  CellType = character(),
  MeanAvgOverlap = numeric(),
  stringsAsFactors = FALSE
)

# Loop through filtered names and compute mean
for (name in tn_names) {
  cell_type <- sub("_t200_rpr_tn", "", name)
  mean_val <- mean(data_list[[name]]$avg_overlap, na.rm = TRUE)
  avg_overlap_summary <- rbind(avg_overlap_summary, data.frame(
    CellType = cell_type,
    MeanAvgOverlap = mean_val
  ))
}

print(avg_overlap_summary)


#####   average ribosomal per cell type  =============================================================
rb_names <- grep("t200_rpr_rb$", names(data_list), value = TRUE)
# Initialize empty list to store results
avg_overlap_summary <- data.frame(
  CellType = character(),
  MeanAvgOverlap = numeric(),
  stringsAsFactors = FALSE
)

# Loop through filtered names and compute mean
for (name in rb_names) {
  cell_type <- sub("_t200_rpr_rb", "", name)
  mean_val <- mean(data_list[[name]]$meanTop200Ctl, na.rm = TRUE)
  avg_overlap_summary <- rbind(avg_overlap_summary, data.frame(
    CellType = cell_type,
    MeanAvgOverlap = mean_val
  ))
}


#####   TF distribution =============================================================
# Get cell types from the names of data_list
cell_types <- sub("_t200_rpr_tf", "", grep("t200_rpr_tf$", names(data_list), value = TRUE))

# Initialize result data frame
tf_vs_tn_summary <- data.frame(
  CellType = character(),
  MaxOverlap = numeric(),
  NumGenes = integer(),
  NumAboveThreshold = integer(),
  PercentAbove = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each cell type and compute the percentage
for (cell in cell_types) {
  # Get corresponding tf and tn datasets
  tf_data <- data_list[[paste0(cell, "_t200_rpr_tf")]]
  tn_data <- data_list[[paste0(cell, "_t200_rpr_tn")]]
  
  # Get the max avg_overlap from tn
  max_overlap <- max(tn_data$avg_overlap, na.rm = TRUE)
  
  # Count how many genes in tf have meanTop200_AD > max_overlap
  num_genes <- nrow(tf_data)
  num_above <- sum(tf_data$meanTop200_AD > max_overlap, na.rm = TRUE)
  
  # Calculate percentage
  percent_above <- (num_above / num_genes) * 100
  
  # Store results
  tf_vs_tn_summary <- rbind(tf_vs_tn_summary, data.frame(
    CellType = cell,
    MaxOverlap = max_overlap,
    NumGenes = num_genes,
    NumAboveThreshold = num_above,
    PercentAbove = percent_above
  ))
}

# View the summary
print(tf_vs_tn_summary)



### empirical p-value calculation =========================================================
# Initialize result list
empirical_pvals_list <- list()

# Loop over all relevant cell types
cell_types <- sub("_t200_rpr_tf", "", grep("t200_rpr_tf$", names(data_list), value = TRUE))

for (cell in cell_types) {
  tf_data <- data_list[[paste0(cell, "_t200_rpr_tf")]]
  tn_data <- data_list[[paste0(cell, "_t200_rpr_tn")]]
  
  null_distribution <- tn_data$avg_overlap
  
  # Compute empirical p-values for each gene
  tf_data$empirical_pval <- sapply(tf_data$meanTop200_AD, function(obs_val) {
    # Count how many null values are >= observed
    (sum(null_distribution >= obs_val, na.rm = TRUE) + 1) / (length(null_distribution) + 1)
  })
  
  empirical_pvals_list[[cell]] <- tf_data
}

summary_df <- data.frame(CellType = character(), PropBelowPval001 = numeric())

for (cell in names(empirical_pvals_list)) {
  df <- empirical_pvals_list[[cell]]
  prop <- mean(df$empirical_pval < 0.001, na.rm = TRUE)
  summary_df <- rbind(summary_df, data.frame(
    CellType = cell,
    PropBelowPval001 = prop
  ))
}






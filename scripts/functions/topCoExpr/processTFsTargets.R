

### 1. from edgeList to edgeList dataframe 
process_file <- function(file_path, TFs, top_n_targets = NULL) {
  
  # 1. Read in edge lists
  df <- fread(file_path)
  colnames(df) <- c("geneA", "geneB", "rank_norm_corr") 
  
  # 2. Filter TFs
  df_filtered <- df %>% filter(geneA %in% names(TFs)) 
  
  # Sanity check: Ensure filtered TFs match expected names
  int <- unique(df_filtered$geneA)[unique(df_filtered$geneA) %in% names(TFs)]
  if (length(unique(df_filtered$geneA)) != length(int)) {
    stop("Sanity check failed: unique geneA count does not match int.")
  }
  
  # Sanity check: Ensure each geneA has the expected number of partners
  freq_table <- as.data.frame(table(df_filtered$geneA))
  if (!all(freq_table$Freq == length(unique(df$geneA)) -1 )) {
    stop("Sanity check failed: Not all geneA have expected co-expressed partners.")
  }
  
  # Extract top N targets if specified
  if (!is.null(top_n_targets)) {
    df_filtered2 <- df_filtered %>%
      group_by(geneA) %>%
      slice_max(order_by = rank_norm_corr, n = top_n_targets, with_ties = FALSE) %>%
      ungroup()
    print("Done processing top targets")
    return(df_filtered2)
  }
  
  print("Done processing all targets")
  return(df_filtered)
}





### 2. from matrix to edgeList dataframe 
extract_top_coexpressed <- function(matrix_path, TFs, top_n_targets = NULL) {
  
  #mat <- readRDS(matrix_path)
  mat <- matrix_path
  
  # Ensure matrix is symmetric
  if (!all(dim(mat) == dim(mat)[1])) {
    stop("The input matrix must be square (gene-by-gene correlation matrix).")
  }
  
  # Ensure TFs are in the matrix
  Tfs <- intersect(TFs, rownames(mat))
  if (length(Tfs) == 0) {
    stop("None of the provided TFs are found in the matrix.")
  }
  
  # Initialize results list
  results <- list()
  
  # Loop through each TF and extract top co-expressed genes
  for (tf in Tfs) {
    # Extract correlations for the given TF
    tf_cor <- mat[tf, ]
    print(tf)
    
    # Check if all values in tf_cor are NA, if so, skip to the next TF
    if (all(is.na(tf_cor) | is.infinite(tf_cor))) {
      print("Moving onto the next TF")
      next
    }
    
    # Remove self-correlation
    tf_cor <- tf_cor[names(tf_cor) != tf]
    
    # Rank by absolute correlation value (descending order)
    ranked_genes <- sort(tf_cor, decreasing = TRUE)
    
    # Select top N targets if specified
    if (!is.null(top_n_targets)) {
      ranked_genes <- head(ranked_genes, top_n_targets)
    }
    
    # Create a dataframe for the TF
    tf_df <- data.frame(
      geneA = tf,
      geneB = names(ranked_genes),
      rank_correlation = ranked_genes,
      row.names = NULL
    )
    
    results[[tf]] <- tf_df
  }
  
  # Combine results into a single dataframe
  final_df <- do.call(rbind, results)
  return(final_df)
}


# 3. from matrix to 1 TF all ranked targets 
rank_gene_partners <- function(coexpr_matrix, gene) {
  # Check if gene is in the matrix
  if (!(gene %in% rownames(coexpr_matrix))) {
    stop("Gene not found in the matrix.")
  }
  
  # Subset the matrix for the given gene
  gene_coexpr <- coexpr_matrix[gene, , drop = FALSE]
  gene_coexpr <- as.data.frame(t(gene_coexpr))
  colnames(gene_coexpr) <- "coExpr_value"
  
  # Remove self-correlation
  gene_coexpr <- gene_coexpr[rownames(gene_coexpr) != gene, , drop = FALSE]
  
  # Rank all partners based on correlation (descending order)
  gene_coexpr$coExpr_rank <- rank(-gene_coexpr$coExpr_value, ties.method = "min")
  
  # Get top 200 partners
  top200_partners <- head(order(-gene_coexpr$coExpr_value), 200)
  
  # Assign top200 rank
  gene_coexpr$top200_rank <- NA
  gene_coexpr$top200_rank[top200_partners] <- rank(-gene_coexpr$coExpr_value[top200_partners], ties.method = "min")
  
  # Indicate whether a partner is in the top 200
  gene_coexpr$top200 <- FALSE
  gene_coexpr$top200[top200_partners] <- TRUE
  
  # Format final dataframe
  result_df <- data.frame(
    geneA = gene,
    geneB = rownames(gene_coexpr),
    coExpr_value = gene_coexpr$coExpr_value, 
    coExpr_rank = gene_coexpr$coExpr_rank,
    top200_rank = gene_coexpr$top200_rank,
    top200 = gene_coexpr$top200
  )
  
  return(result_df)
}

#rank_gene_partners (readRDS("/home/amorin/scratch/aggregate_microglia/Cormats/Hg_pcor/aggregate_cormat_allrank_filter_lenient_hg.RDS")[["Agg_mat"]],"RUNX1")


# 4. TF-gene average rank coexpression =================================================
average_tf_target_ranks <- function(mic_TF_allTargets, condition, tf_by_gene_binding) {

  # Step 1: Get relevant TF-target sets by condition
  suffix <- paste0("_", condition)
  filtered_keys <- names(mic_TF_allTargets)[endsWith(names(mic_TF_allTargets), suffix)]
  
  # Define sets for filtering
  valid_geneA <- intersect(colnames(tf_by_gene_binding), colnames(tf_by_gene_binding))
  valid_geneB <- tf_by_gene_binding$geneB

  cat("Valid geneA count (TFs):", length(valid_geneA), "\n")
  cat("Valid geneB count (targets):", length(valid_geneB), "\n")
  
  # Step 2: Filter each list element
  print("filtering datasets elements")
  filtered_list <- lapply(filtered_keys, function(k) {
    df <- mic_TF_allTargets[[k]]
    df %>%
      filter(geneA %in% valid_geneA, geneB %in% valid_geneB)
  })
  
  # Combine all into one dataframe
  print("combining all datasets elements")
  combined_df <- bind_rows(filtered_list)
  
  # Step 3: Average rank_norm_corr per geneA-geneB pair
  print("averaging TF-gene coexpression ranking")
  summary_df <- combined_df %>%
    group_by(geneA, geneB) %>%
    summarize(
      avg_rank_norm_corr = mean(rank_norm_corr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(condition = condition)
  
  # Add binding score
  cat("🔗 Adding binding scores...\n")
  
  # Reshape to long format
  tf_binding_long <- tf_by_gene_binding %>%
    pivot_longer(-geneB, names_to = "geneA", values_to = "binding_score")
  
  summary_df <- summary_df %>%
    left_join(tf_binding_long, by = c("geneA", "geneB"))
  
  cat("🎉 TF-target rank summarization completed successfully with binding scores included!\n")
  return(summary_df)

}


# 5. TF-gene aggregate coexpression =================================================


tf_gene_aggCoExpr_litCur <- function(mic_TF_allTargets, condition, lit_targets) {
  cat("starting TF-target aggregation for condition:", condition, "\n")

  suffix <- paste0("_", condition)
  dataset_keys <- names(mic_TF_allTargets)[endsWith(names(mic_TF_allTargets), suffix)]
  cat("Found", length(dataset_keys), "datasets matching suffix", suffix, "\n")

  if (length(dataset_keys) == 0) {
    stop("❌ No datasets match the given condition.")
  }

  tf_set <- unique(lit_targets$TF_Symbol)
  target_set <- unique(lit_targets$Target_Symbol)

  cat("Literature-curated TFs:", length(tf_set), "\n")
  cat("Literature-curated Targets:", length(target_set), "\n")

  # Filter each list element
  cat("filtering dataset elements") 
  filtered_list <- lapply(dataset_keys, function(k) {
    df <- mic_TF_allTargets[[k]]
    initial_n <- nrow(df)
    df_filtered <- df %>%
      filter(geneA %in% tf_set, geneB %in% target_set)
    #cat("Dataset:", k, "→", nrow(df_filtered), "rows retained from", initial_n, "\n")
    return(df_filtered)
  })

  combined_df <- bind_rows(filtered_list)
  cat("Total combined TF-target pairs:", nrow(combined_df), "\n")

  if (nrow(combined_df) == 0) {
    stop("❌ No geneA/geneB pairs matched literature-curated targets.")
  }

  # Aggregate rank_norm_corr across all datasets
  summary_df <- combined_df %>%
    group_by(geneA, geneB) %>%
    summarize(
      avg_rank_norm_corr = mean(rank_norm_corr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(condition = condition)

  cat("Aggregated TF-target ranks:", nrow(summary_df), "\n")

  # Determine if TF-target pair is curated
  lit_pairs <- lit_targets %>%
    select(TF_Symbol, Target_Symbol) %>%
    distinct() %>%
    mutate(key = paste(TF_Symbol, Target_Symbol, sep = "_"))

  summary_df <- summary_df %>%
    mutate(key = paste(geneA, geneB, sep = "_")) %>%
    mutate(TF_target_predicted = if_else(key %in% lit_pairs$key, "positive", "negative")) %>%
    select(geneA, geneB, rank_norm_corr, condition, TF_target_predicted)

  cat("🎉 TF-target summary completed. Positive predictions:",
      sum(summary_df$TF_target_predicted == "positive"), "\n")

  return(summary_df)
}

# 6. add other mddalities =================================================
# A. aggregate TF-gene coexpression profiles for all pairs 
# B. integrate binding scores for the avaibale gene pairs 
# C. integrate literature curation evidence for available gene pairs 

library(tidyr)
library(tidyverse)

summarize_TF_target_profiles <- function(mic_TF_allTargets, condition, tf_by_gene_binding, lit_targets) {
  cat("starting TF-target summarization for condition:", condition, "\n")
  
  # Step 1: Filter datasets by condition suffix
  suffix <- paste0("_", condition)
  relevant_keys <- names(mic_TF_allTargets)[endsWith(names(mic_TF_allTargets), suffix)]
  cat("Found", length(relevant_keys), "datasets with suffix", suffix, "\n")
  
  if (length(relevant_keys) == 0) stop(" No datasets match the given condition.")
  
  # Step 2: Combine all relevant datasets and tag dataset ID
  combined_df <- bind_rows(lapply(relevant_keys, function(k) {
    df <- mic_TF_allTargets[[k]]
    df$dataset_id <- k
    return(df)
  }))
  
  cat(" Combined", nrow(combined_df), "rows from", length(relevant_keys), "datasets\n")
  
  # Step 3a: Keep per-dataset rank_norm_corr
  per_dataset_df <- combined_df %>%
    select(geneA, geneB, dataset_id, rank_norm_corr) %>%
    mutate(condition = condition)
  
  # Step 3b: Aggregate across datasets
  summary_df <- combined_df %>%
    group_by(geneA, geneB) %>%
    summarize(avg_rank_norm_corr = mean(rank_norm_corr, na.rm = TRUE), .groups = "drop") %>%
    mutate(condition = condition)
  
  cat("Aggregated coexpression for", nrow(summary_df), "unique gene pairs\n")
  
  # Step 4: Add binding scores
  cat("Adding binding scores...\n")
  tf_binding_long <- tf_by_gene_binding %>%
    pivot_longer(-geneB, names_to = "geneA", values_to = "binding_score")
  
  summary_df <- summary_df %>%
    left_join(tf_binding_long, by = c("geneA", "geneB"))
  
  # Step 5: Add literature support
  cat("Mapping literature-supported pairs...\n")
  lit_pairs <- lit_targets %>%
    distinct(TF_Symbol, Target_Symbol) %>%
    mutate(pair_key = paste(TF_Symbol, Target_Symbol, sep = "_"))
  
  summary_df <- summary_df %>%
    mutate(pair_key = paste(geneA, geneB, sep = "_")) %>%
    mutate(
      TF_target_predicted = if_else(pair_key %in% lit_pairs$pair_key, "positive", "negative")
    )
  
  # Step 6: Attach per-dataset values back into summary
  final_df <- summary_df %>%
    left_join(per_dataset_df, by = c("geneA", "geneB", "condition"))
  
  # Order columns nicely
  final_df <- final_df %>%
    select(geneA, geneB, dataset_id, rank_norm_corr, avg_rank_norm_corr,
           condition, binding_score, TF_target_predicted)
  
  cat("Summary complete. Literature-positive pairs:",
      sum(final_df$TF_target_predicted == "positive", na.rm = TRUE), "\n")
  
  return(final_df)
}




summarize_TF_target_profilesV1 <- function(tf_top200_list, condition, cell_type,
                                         tf_by_gene_binding, lit_targets, tf_rpr) {
  cat("Starting TF-target summarization for:", cell_type, "-", condition, "\n")
  
  # Step 1: Extract relevant datasets for the given cell type
  if (!cell_type %in% names(tf_top200_list)) {
    stop("Provided cell type not found in tf_top200_list.")
  }
  
  celltype_list <- tf_top200_list[[cell_type]]
  
  # Filter datasets by condition suffix (e.g., "_AD", "_CTL")
  suffix <- paste0("_", condition)
  relevant_keys <- names(celltype_list)[endsWith(names(celltype_list), suffix)]
  cat("Found", length(relevant_keys), "datasets for cell type", cell_type, "with condition", condition, "\n")
  
  if (length(relevant_keys) == 0) stop("No datasets found for the specified cell type and condition.")
  
  # Step 2: Combine relevant data
  combined_df <- bind_rows(lapply(relevant_keys, function(k) {
    df <- celltype_list[[k]]
    df$dataset_id <- k
    return(df)
  }))
  
  cat("Combined", nrow(combined_df), "rows from", length(relevant_keys), "datasets\n")
  
  # Step 3: Average correlation ranks
summary_df <- combined_df %>%
  group_by(geneA, geneB) %>%
  summarize(
    avg_rank_norm_corr = mean(rank_norm_corr, na.rm = TRUE),
    n_datasets = n_distinct(dataset_id),  # Count of distinct datasets
    .groups = "drop"
  ) %>%
  mutate(
    condition = condition,
    cell_type = cell_type
  )
  
  cat("Aggregated coexpression for", nrow(summary_df), "unique gene pairs\n")
  
  # Step 4: Add binding scores (long format)
  cat("Adding binding scores...\n")
  tf_binding_long <- tf_by_gene_binding %>%
    pivot_longer(cols = -geneB, names_to = "geneA", values_to = "binding_score")
  
  summary_df <- summary_df %>%
    left_join(tf_binding_long, by = c("geneA", "geneB"))
  
  # Step 5: Add literature support
  cat("Mapping literature-supported pairs...\n")
  lit_pairs <- lit_targets %>%
    distinct(TF_Symbol, Target_Symbol) %>%
    mutate(pair_key = paste(TF_Symbol, Target_Symbol, sep = "_"))
  
  summary_df <- summary_df %>%
    mutate(pair_key = paste(geneA, geneB, sep = "_")) %>%
    mutate(
      TF_target_predicted = if_else(pair_key %in% lit_pairs$pair_key, "positive", "negative")
    )
  
  # Step 6: Add TF reproducibility scores
  cat("Adding TF reproducibility scores...\n")
  tf_rpr <- tf_rpr %>%
    group_by(cell_type) %>%
    mutate(
      mean_ranked_rpr = rank(meanTop200Ctl, ties.method = "average"),
      mean_ranked_rpr_scaled = (mean_ranked_rpr - min(mean_ranked_rpr)) / 
        (max(mean_ranked_rpr) - min(mean_ranked_rpr))
    ) %>%
    ungroup()

  tf_rpr_clean <- tf_rpr %>%
    filter(cell_type == !!cell_type) %>%
    select(geneA, mean_ranked_rpr_scaled) %>%
    rename(tf_rpr_score = mean_ranked_rpr_scaled)
  
  summary_df <- summary_df %>%
    left_join(tf_rpr_clean, by = "geneA")
  
  # Step 7: Final formatting
  final_df <- summary_df %>%
    select(geneA, geneB, avg_rank_norm_corr, n_datasets, condition, cell_type, binding_score,
           TF_target_predicted, tf_rpr_score)
  
  cat("Summary complete. Literature-positive pairs:",
      sum(final_df$TF_target_predicted == "positive", na.rm = TRUE), "\n")
  
  return(final_df)
}

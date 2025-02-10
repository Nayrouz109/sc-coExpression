

### from edgeList to edgeList dataframe 
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





### from matrix to edgeList dataframe 
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







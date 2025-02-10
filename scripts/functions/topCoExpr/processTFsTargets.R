

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












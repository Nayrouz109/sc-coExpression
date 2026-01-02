

#### these are all functions used for the brain cell type coexpression conservation (mouse) analysis 

library(tidyr)

### 1. convert mouse names into humans 
convert_mouse_to_human <- function(mouse_cor_mat, mouse_to_human) {
    
    # Ensure that mouse_to_human is a named vector
    if (is.null(names(mouse_to_human))) {
        stop("mouse_to_human must be a named vector with mouse gene symbols as names")
    }

    # Convert row and column names to human symbols
    human_genes <- mouse_to_human[rownames(mouse_cor_mat)]
  
    # Remove genes with NA human orthologues
    valid_genes <- !is.na(human_genes)
    filtered_mat <- mouse_cor_mat[valid_genes, valid_genes]
  
    # Rename the matrix with human gene symbols
    rownames(filtered_mat) <- human_genes[valid_genes]
    colnames(filtered_mat) <- human_genes[valid_genes]
  
    return(filtered_mat)
}
# Example usage
# filtered_matrix <- convert_mouse_to_human(alex_mouse_all_rank$Agg_mat, mouse_to_human)



### 2. extract all partners per genes  
library(data.table)
library(dplyr)

process_file <- function(df, TFs = NULL, top_n_targets = NULL) {
  
  # 1. Read in edge lists
  #df <- fread(file_path)
  colnames(df) <- c("geneA", "geneB", "rank_norm_corr") 
  
  # 2. Filter TFs if provided
  if (!is.null(TFs)) {
    df_filtered <- df %>% filter(geneA %in% TFs) 
    # Sanity check: Ensure all filtered TFs exist in the provided list
    if (!all(unique(df_filtered$geneA) %in% TFs)) {
        stop("Sanity check failed: Some geneA values are not in the provided TF list.")
    }
  } else {
    df_filtered <- df
  }
  
  # 3. Sanity check: Ensure each geneA has at least one co-expressed partner
  freq_table <- as.data.frame(table(df_filtered$geneA))
  if (any(freq_table$n == 0)) {
    stop("Sanity check failed: Some geneA have no co-expressed partners.")
  }
  
  # 4. Extract top N targets if specified
  if (!is.null(top_n_targets)) {
    df_filtered2 <- df_filtered %>%
      group_by(geneA) %>%
      slice_max(order_by = rank_norm_corr, n = top_n_targets, with_ties = FALSE) %>%
      ungroup()
    message("Done processing top targets")
    return(df_filtered2)
  }
  
  message("Done processing all targets")
  return(df_filtered)
}



### this will give you the top 200 parterns per TF across all studies (the most frequent partners) 
library(purrr)
process_cellType_top200 <- function(cellType_top200, avg = NULL) {
  cellType_top200_summary <- bind_rows(
    lapply(names(cellType_top200), function(tf) {
      df <- cellType_top200[[tf]]
      
      if (is.data.frame(df)) {
        # If avg is TRUE, select top 200 geneB based on freq_CTL
        if (!is.null(avg) && avg == TRUE) {
          df <- df %>% arrange(desc(freq_CTL)) %>% slice_head(n = 200)
        }
        
        df %>%
          select(geneB) %>%
          mutate(TF = tf)
      } else {
        NULL  
      }
    })
  ) %>% select(TF, geneB)

  return(cellType_top200_summary)
}
# Example usage
# mic_top200_summary <- process_cellType_top200(cellType_top200, avg = TRUE)




# this is to filter out human edge list and mouse cor matrix to match gene counts 
filter_edge_lists_and_matrix <- function(input_files, mouse_matrix) {
    library(data.table)
    # Read edge lists and store them in a named list
    print("reading files")
    edge_lists <- lapply(input_files, fread)
    names(edge_lists) <- basename(input_files)  # Assign file names as list names
  
    # Extract unique genes from edge lists
    edge_genes <- unique(unlist(lapply(edge_lists, function(df) unique(c(df$source, df$target)))))
  
    # Extract gene names from the mouse correlation matrix
    matrix_genes <- rownames(mouse_matrix)
  
    # Find common genes
    common_genes <- intersect(edge_genes, matrix_genes)
  
    # Filter edge lists to retain only interactions with common genes
    print("filtering hg mtx")
    filtered_edge_lists <- lapply(edge_lists, function(df) {
        df[df$source %in% common_genes & df$target %in% common_genes, ]
    })
  
    # Filter mouse correlation matrix to retain only common genes
    print("filtering mm mtx")
    filtered_matrix <- mouse_matrix[common_genes, common_genes, drop = FALSE]
  
    return(list(filtered_edges = filtered_edge_lists, filtered_matrix = filtered_matrix))
}

# Example usage:
#filtered_data <- filter_edge_lists_and_matrix(input_files, flt_mm_all_rank)
# Access filtered edge lists and matrix
#filtered_edges <- filtered_data$filtered_edges
#filtered_matrix <- filtered_data$filtered_matrix



### 4. generate gene-target summaries 
integrate_gene_pairs_summaries <- function(gn_allTrg_list, TFs, freqTable, condition, ribo_list, matched_data_mm) {
    
    library(dplyr)
    
    # Filter studies based on the condition
    relevant_studies <- names(gn_allTrg_list)[grepl(condition, names(gn_allTrg_list))]
  
    # Combine relevant studies into one dataframe
    combined_df <- bind_rows(gn_allTrg_list[relevant_studies], .id = "study")
  
    # Compute avg_cor for each geneA-geneB pair
    print("computing hg avg cor")
    avg_cor_df <- combined_df %>%
        group_by(geneA, geneB) %>%
        summarize(hg_avg_cor = mean(rank_norm_corr, na.rm = TRUE), .groups = "drop")
  
    # Determine if geneA is a TF
    avg_cor_df <- avg_cor_df %>%
        mutate(is_TF = geneA %in% TFs) %>%
        mutate(is_ribo = geneA %in% ribo_list & geneB %in% ribo_list)
    
    # Determine if geneB is among the top 200 in freqTable for geneA
    top200_table <- freqTable %>%
        arrange(desc(freq_CTL)) %>%  # Sort by freq_CTL descending
        group_by(geneA) %>%
        slice_head(n = 200) %>%  # Take the top 200 per geneA
        ungroup() %>%
        select(geneA, geneB) %>%
        mutate(is_top200_hg = TRUE)
  
    # Merge with avg_cor_df to annotate top 200 geneB
    print("merging top200 info")
    final_table <- avg_cor_df %>%
        left_join(top200_table, by = c("geneA", "geneB")) %>%
        mutate(is_top200_hg = replace_na(is_top200_hg, FALSE))
    
    print("extracting mouse cor")
    matched_data_mm <- matched_data_mm %>%
        rename(mm_avg_cor = cor, is_top200_mm = is_top200)
    
    final_table <- final_table %>%
        left_join(matched_data_mm, by = c("geneA", "geneB"))
  
    return(final_table)
}

# Example usage:
#summary_table <- integrate_gene_pairs_summaries(gn_allTrg_list, TFs, freqTable, "AD", ribo_list)



### from matrix to edge list is_top200
generate_edge_list <- function(cor_matrix, top_n_targets = 200) {
  # Convert correlation matrix into a long-format data frame (edge list)
  edge_list <- as.data.frame(as.table(cor_matrix))
  colnames(edge_list) <- c("geneA", "geneB", "cor")
  
  # Ensure geneA and geneB are character columns
  edge_list$geneA <- as.character(edge_list$geneA)
  edge_list$geneB <- as.character(edge_list$geneB)
  
  # Filter out self-correlations
  edge_list <- edge_list %>%
    filter(geneA != geneB)
  
  # Identify the top 200 correlated partners for each geneA
  top200_genes <- edge_list %>%
    group_by(geneA) %>%
    slice_max(order_by = cor, n = top_n_targets, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(is_top200 = TRUE) %>%
    select(geneA, geneB, is_top200)
  
  # Merge with the original edge list
  edge_list <- edge_list %>%
    left_join(top200_genes, by = c("geneA", "geneB")) %>%
    mutate(is_top200 = replace_na(is_top200, FALSE))
  
  return(edge_list)
}



compute_geneA_summary <- function(intg_summary) {
  library(dplyr)
  intg_summary %>%
    group_by(geneA) %>%
    summarize(
      s.cor = cor(hg_avg_cor, mm_avg_cor, method = "spearman"),
      top_200 = sum(is_top200_hg == TRUE & is_top200_mm == TRUE)
    )
}

compute_geneA_summary <- function(intg_summary) {
  intg_summary %>%
    group_by(geneA) %>%
    summarize(
      s.cor = ifelse(sd(hg_avg_cor) == 0 | sd(mm_avg_cor) == 0, NA, cor(hg_avg_cor, mm_avg_cor, method = "spearman")),
      top_200 = sum(is_top200_hg == TRUE & is_top200_mm == TRUE)
    )
}

# Example usage:
# result <- compute_geneA_summary(intg_summary)
# head(result)



library(data.table)
library(dplyr)
library(purrr)
library(magrittr)





# ------------------ compute cor differences ---------------------------------------

signed_rank_normalize <- function(x) {
  # Rank normalize separately for positive and negative values
  pos_ranks <- rank(x[x >= 0], ties.method = "min") / length(x)
  neg_ranks <- -rank(-x[x < 0], ties.method = "min") / length(x)
  
  # Combine positive and negative ranks while maintaining original order
  normalized_x <- numeric(length(x))
  normalized_x[x >= 0] <- pos_ranks
  normalized_x[x < 0] <- neg_ranks
  
  return(normalized_x)
}

compute_signed_rank_norm_diff <- function(mic_list, ast_list, exc_list) {
  # Extract only control datasets
  mic_ctl <- mic_list[grep("_CTL$", names(mic_list))]
  ast_ctl <- ast_list[grep("_CTL$", names(ast_list))]
  exc_ctl <- exc_list[grep("_CTL$", names(exc_list))]
  
  study_names <- intersect(names(mic_ctl), intersect(names(ast_ctl), names(exc_ctl)))
  
  result_list <- lapply(study_names, function(study) {
    
    print(paste("working on" , study))
    mic_df <- mic_ctl[[study]]
    ast_df <- ast_ctl[[study]]
    exc_df <- exc_ctl[[study]]
    
    # Merge by gene pair (geneA, geneB)
    merged_ast <- merge(mic_df, ast_df, by = c("geneA", "geneB"), suffixes = c("_mic", "_ast"))
    merged_exc <- merge(mic_df, exc_df, by = c("geneA", "geneB"), suffixes = c("_mic", "_exc"))
    merged_ea  <- merge(ast_df, exc_df, by = c("geneA", "geneB"), suffixes = c("_ast", "_exc"))
    
    # Compute signed differences in rank_norm_corr
    print(paste("computing cor diff"))
    merged_ast$signed_diff <- merged_ast$rank_norm_corr_mic - merged_ast$rank_norm_corr_ast
    merged_exc$signed_diff <- merged_exc$rank_norm_corr_mic - merged_exc$rank_norm_corr_exc
    merged_ea$signed_diff  <- merged_ea$rank_norm_corr_ast - merged_ea$rank_norm_corr_exc
    
    # Rank normalize differences while preserving the sign
    print(paste("ranking cor diff"))
    merged_ast$rank_norm_diff <- signed_rank_normalize(merged_ast$signed_diff)
    merged_exc$rank_norm_diff <- signed_rank_normalize(merged_exc$signed_diff)
    merged_ea$rank_norm_diff  <- signed_rank_normalize(merged_ea$signed_diff)
    
    # Return a dataframe with relevant columns
    data.frame(
      geneA = merged_ast$geneA,
      geneB = merged_ast$geneB,
      rank_norm_diff_ast = merged_ast$rank_norm_diff,
      rank_norm_diff_exc = merged_exc$rank_norm_diff,
      rank_norm_diff_ea  = merged_ea$rank_norm_diff ,
      study = study
    )
  })
  
  # Combine all study results into one dataframe
  print(paste("combining study result"))
  combined_df <- do.call(rbind, result_list)
  
  # Compute the average rank-normalized difference across studies per gene pair
  print(paste("aggregating data"))
  summary_df <- aggregate(. ~ geneA + geneB, data = combined_df[, c("geneA", "geneB", "rank_norm_diff_ast", "rank_norm_diff_exc", "rank_norm_diff_ea")], mean)
  
  return(summary_df)
}





# Function to extract rank_norm_corr values
extract_rank_norm_corr_ribo <- function(cell_type_lists, geneList) {
  results_list <- list()
  
  for (cell_type in names(cell_type_lists)) {
    print(paste("Processing cell type:", cell_type))
    df_list <- cell_type_lists[[cell_type]] %>% mutate( 
                          CellType = cell_type, 
                          #gene_pair = paste(geneA, geneB, sep = "_"),
                          dataset = sub("^[^_]+_", "", dataset), 
                          gene_pair_ds = paste0(dataset, "_" ,gene_pair)
    ) %>% rowwise() %>% mutate(
                gene_pair = paste(sort(c(geneA, geneB)), collapse = "_")
    ) %>%
    ungroup()

    

    #print("extracting cor values ... ")
    extracted_values <- df_list %>% filter(geneA %in% geneList & geneB %in% geneList) 
    #print(head(extracted_values))
    
    if (nrow(extracted_values) > 0) {
      results_list[[cell_type]] <- extracted_values
    }
  }
  
  final_result <- do.call(rbind, results_list)
  print(head(final_result))
  return(final_result)
}


extract_rank_norm_corrV1 <- function(cell_type_lists, top_bottom_combined) {
  results_list <- list()
  
  for (cell_type in names(cell_type_lists)) {
    print(paste("Processing cell type:", cell_type))
    df_list <- cell_type_lists[[cell_type]] %>% mutate( 
                          CellType = cell_type, 
                          gene_pair = paste(geneA, geneB, sep = "_"),
                          dataset = sub("^[^_]+_", "", dataset), 
                          gene_pair_ds = paste0(dataset, "_" ,gene_pair)
    )
    

    print("extracting cor values ... ")
    extracted_values <- df_list %>% filter(gene_pair %in% top_bottom_combined$gene_pair) 
    print("extracting cor category ... ")
    extracted_values <- left_join(extracted_values, top_bottom_combined %>% dplyr::select(gene_pair, category), by = "gene_pair")
    
    if (nrow(extracted_values) > 0) {
      results_list[[cell_type]] <- extracted_values
    }
  }
  
  final_result <- do.call(rbind, results_list)
  return(final_result)
}




# Function to extract rank_norm_corr values

extract_rank_norm_corr <- function(cell_type_lists, geneAA, geneBB) {
  results_list <- list()
  
  for (cell_type in names(cell_type_lists)) {
    df_list <- cell_type_lists[[cell_type]]
    
    print(paste("Processing cell type:", cell_type))
    
    extracted_values <- list()
    for (i in 1:length(geneAA)) {
      filtered_data <- df_list %>% filter((geneA == geneAA[i] & geneB == geneBB[i]) 
                                          #| #(geneA == geneBB[i] & geneB == geneAA[i])
      ) %>% mutate(CellType = cell_type) 
      extracted_values[[i]] <- filtered_data
    }
    extracted_values <- do.call(rbind, extracted_values)
    
    #extracted_values <- df_list[(df_list$geneA %in% geneAA & df_list$geneB %in% geneBB) | 
    #                    (df_list$geneA %in% geneBB & df_list$geneB %in% geneAA), ] %>% mutate(CellType = cell_type) 
    
    
    if (!is.null(extracted_values)) {
      results_list[[cell_type]] <- extracted_values
    }
  }
  
  final_result <- do.call(rbind, results_list)
  return(final_result)
}

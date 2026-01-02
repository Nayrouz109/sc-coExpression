


##### freq tables  ===================================================================================

# 0. edge list of top200 to freq table across conditions across studies
summarize_gene_pairs <- function(edge_lists) {
  library(dplyr)
  
  # 1. Identify unique gene pairs across all studies
  all_pairs <- bind_rows(edge_lists, .id = "study") %>%
    select(geneA, geneB) %>%
    distinct()
  
  # 2. Initialize result dataframe
  summary_df <- all_pairs %>%
    mutate(freq_CTL = 0, freq_AD = 0)
  
  # 3. Iterate through studies and count occurrences
  for (study in names(edge_lists)) {
    df <- edge_lists[[study]]
    condition <- ifelse(grepl("_CTL$", study), "CTL", "AD")
    #print(study)
    #print(condition)
    counts <- df %>%
      dplyr::count(geneA, geneB, name = "freq")
    
    # Merge counts into summary_df
    if (condition == "CTL") {
      summary_df <- summary_df %>%
        left_join(counts, by = c("geneA", "geneB")) %>%
        mutate(freq_CTL = ifelse(is.na(freq), freq_CTL, freq_CTL + freq)) %>%
        select(-freq)
    } else {
      summary_df <- summary_df %>%
        left_join(counts, by = c("geneA", "geneB")) %>%
        mutate(freq_AD = ifelse(is.na(freq), freq_AD, freq_AD + freq)) %>%
        select(-freq)
    }
  }
  
  summary_df <- summary_df %>%
    mutate(delta_freq = freq_AD - freq_CTL) %>% 
    mutate(delta_abs = abs(delta_freq))
  
  return(summary_df)
}

# Example usage:
# result_df <- summarize_gene_pairs(my_data_TF_top200)



summarize_tf_targetsV0 <- function(summary_df) {
  library(dplyr)
  
  tf_summary <- summary_df %>%
    group_by(geneA) %>%
    summarize(
      num_reproducible_targetsCTL = sum(freq_CTL == 5),
      num_reproducible_targetsAD = sum(freq_AD == 5),
      total_targets = n(),
      proportion_reproducible_CTL = num_reproducible_targetsCTL / total_targets,
      proportion_reproducible_AD = num_reproducible_targetsAD / total_targets,
      delta_proportion_reproducible = proportion_reproducible_AD - proportion_reproducible_CTL
    ) %>%
    mutate(
      rank_delta_proportion = rank(-delta_proportion_reproducible, ties.method = "min"),
      standardized_rank = (rank_delta_proportion - min(rank_delta_proportion)) / 
        (max(rank_delta_proportion) - min(rank_delta_proportion))
    ) %>%
    arrange(desc(delta_proportion_reproducible))
  
  return(tf_summary)
}



summarize_tf_targets <- function(summary_df) {
  library(dplyr)
  
  tf_summary <- summary_df %>%
    group_by(geneA) %>%
    summarize(
      num_reproducible_targetsCTL = sum(freq_CTL >= 5),
      num_reproducible_targetsAD = sum(freq_AD >= 5),
      total_targets_CTL = sum(freq_CTL > 0),
      total_targets_AD = sum(freq_AD > 0),
      proportion_reproducible_CTL = num_reproducible_targetsCTL / total_targets_CTL,
      proportion_reproducible_AD = num_reproducible_targetsAD / total_targets_AD,
      delta_proportion_reproducible = proportion_reproducible_CTL - proportion_reproducible_AD 
    ) %>%
    mutate(
      rank_delta_proportion = rank(-delta_proportion_reproducible, ties.method = "min"),
      standardized_rank = (rank_delta_proportion - min(rank_delta_proportion)) / 
        (max(rank_delta_proportion) - min(rank_delta_proportion))
    ) %>%
    arrange(desc(delta_proportion_reproducible))
  
  return(tf_summary)
}

# tf_summary_df <- summarize_tf_targets(result_df)

# ---- TF reproducability ----
##### TF reproducability ============================================================================================

calculate_tf_target_reproducibility <- function(data_list) {

  ad_dfs <- data_list[grep("_AD$", names(data_list))]
  print(length(ad_dfs))
  ctl_dfs <- data_list[grep("_CTL$", names(data_list))]
  print(length(ctl_dfs))
  
  # Function to compute top 200 overlap for a given geneA across datasets
  compute_avg_overlap <- function(df_list) {
    geneA_list <- unique(unlist(lapply(df_list, function(df) df$geneA)))
    
    avg_overlap <- map_dfr(geneA_list, function(gene) {
      geneB_sets <- lapply(df_list, function(df) {
        df %>% filter(geneA == gene) %>% top_n(200, rank_norm_corr) %>% pull(geneB)
      })
      
      # Compute pairwise overlaps
      if (length(geneB_sets) > 1) {
        overlap_values <- combn(geneB_sets, 2, function(pair) {
          length(intersect(pair[[1]], pair[[2]])) #/ 200
        })
        avg_overlap <- mean(overlap_values, na.rm = TRUE)
      } else {
        avg_overlap <- NA
      }
      
      tibble(geneA = gene, meanTop200 = avg_overlap)
    })
    
    return(avg_overlap)
  }
  
  # Compute for AD and Control
  ad_results <- compute_avg_overlap(ad_dfs) %>% rename(meanTop200_AD = meanTop200)
  ctl_results <- compute_avg_overlap(ctl_dfs) %>% rename(meanTop200Ctl = meanTop200)
  
  # Merge results
  final_result <- full_join(ad_results, ctl_results, by = "geneA")
  
  return(final_result)
}


compute_null_overlap <- function(data_list, n) {
  # Extract control datasets
  ctl_dfs <- data_list[grep("_CTL$", names(data_list))]
  
  # Function to compute pairwise overlap
  compute_overlap <- function(selected_targets) {
    overlap_matrix <- sapply(selected_targets, function(x) {
      sapply(selected_targets, function(y) {
        #length(intersect(x, y)) / length(union(x, y))
        length(intersect(x, y)) #/ length(union(x, y))
      })
    })
    mean(overlap_matrix[lower.tri(overlap_matrix)]) # Average pairwise overlap
  }
  
  # Initialize result storage
  results <- numeric(n)
  
  # Perform n iterations
  for (i in 1:n) {
    
    print(paste("Iteration:", i)) # Print iteration number
    
    # Select 1 random TF per dataset and extract its top 200 targets
    selected_targets <- lapply(ctl_dfs, function(df) {
      sampled_tf <- sample(unique(df$geneA), 1) # Pick a random TF
      df[df$geneA == sampled_tf, ]$geneB # Extract top 200 targets
    })
    
    # Compute average pairwise overlap
    results[i] <- compute_overlap(selected_targets)
  }
  
  # Create a summary table
  summary_table <- data.frame(iteration = 1:n, avg_overlap = results)
  
  return(summary_table)
}


##### TF reproducability (Leave-One-Study-Out Analysis) ============================================================================================
# 1. Recompute the reproducibility scores by removing one study at a time.
# 2. Compare each "leave-one-out" result to the full dataset.
# 3. Flag studies that significantly alter the results when removed.

calculate_tf_target_reproducibility_loso <- function(data_list) {
  
  # Identify unique study names without _AD or _CTL suffix
  study_names <- unique(gsub("_(AD|CTL)$", "", names(data_list)))
  
  # Function to compute top 200 overlap for a given geneA across datasets
  compute_avg_overlap <- function(df_list) {
    geneA_list <- unique(unlist(lapply(df_list, function(df) df$geneA)))
    
    avg_overlap <- map_dfr(geneA_list, function(gene) {
      geneB_sets <- lapply(df_list, function(df) {
        df %>% filter(geneA == gene) %>% top_n(200, rank_norm_corr) %>% pull(geneB)
      })
      
      if (length(geneB_sets) > 1) {
        overlap_values <- combn(geneB_sets, 2, function(pair) {
          length(intersect(pair[[1]], pair[[2]])) 
        })
        avg_overlap <- mean(overlap_values, na.rm = TRUE)
      } else {
        avg_overlap <- NA
      }
      
      tibble(geneA = gene, meanTop200 = avg_overlap)
    })
    
    return(avg_overlap)
  }
  
  # Compute full dataset results
  print("computing full ad")
  full_ad_results <- compute_avg_overlap(data_list[grep("_AD$", names(data_list))]) %>%
    rename(meanTop200_AD = meanTop200)
  print("computing full ctl")
  full_ctl_results <- compute_avg_overlap(data_list[grep("_CTL$", names(data_list))]) %>%
    rename(meanTop200_CTL = meanTop200)
  
  # Perform Leave-One-Study-Out analysis
  loso_results <- map_dfr(study_names, function(study) {
    # Identify datasets to keep
    cat(study, "was removed.\n")
    loso_list <- data_list[!grepl(paste0("^", study, "_"), names(data_list))]
    
    ad_reduced <- loso_list[grep("_AD$", names(loso_list))]
    ctl_reduced <- loso_list[grep("_CTL$", names(loso_list))]
    
    ad_reduced_results <- compute_avg_overlap(ad_reduced) %>%
      rename(meanTop200_AD_LOSO = meanTop200)
    
    ctl_reduced_results <- compute_avg_overlap(ctl_reduced) %>%
      rename(meanTop200_CTL_LOSO = meanTop200)
    
    combined_results <- full_join(ad_reduced_results, ctl_reduced_results, by = "geneA") %>%
      mutate(Study_Removed = study)
    
    return(combined_results)
  })
  
  # Merge full results with LOSO results for comparison
  final_result <- full_join(full_ad_results, full_ctl_results, by = "geneA") %>%
    left_join(loso_results, by = "geneA")
  
  return(final_result)
}










##### freq tables + extra info integration  ===================================================================================

### this is to generate integrated ranking dataframe for a given TF, function takes in:
# 1- freqTable                 # my analysis
# 2- tf_targets_alex_fishZ     # alex's mic FZ
# 3- tf_targets_alex_allRank   # alex's mic allRank
# 4- TFs                       # alex's global 
# 5- gene (example "RUNX1")    # gene of interest 
# 6- lit_TFs                   # literature curated targets 
# 7- delta threshold           # could be a number of none 
# 8- mic_neuroinflammation     # neuroinflammation paper
# 9- SEA-AD_mic                # genes mentioned in SEA-AD paper related to MIC
# 10- eRregulon                # GRN from the SEA-AD paper 


integrate_gene_ranking <- function(freqTable, 
                                   tf_targets_alex_fishZ, 
                                   tf_targets_alex_allRank, 
                                   TFs, 
                                   gene, 
                                   lit_TFs, 
                                   delta_threshold = NULL, 
                                   mic_neuroinflammation = NULL, 
                                   SEA_AD_micGenes = NULL, 
                                   eRegulon = NULL) {
  
  # Filter freqTable for the gene of interest
  gene_merged <- freqTable %>% filter(geneA == gene) %>%
    left_join(
      tf_targets_alex_fishZ %>%
        select(geneB, coExpr_rank) %>%
        rename_with(~ paste0(., "_FZ"), c("coExpr_rank")),
      by = "geneB"
    ) %>%
    left_join(
      tf_targets_alex_allRank %>%
        select(geneB, coExpr_rank) %>%
        rename_with(~ paste0(., "_allRank"), c("coExpr_rank")),
      by = "geneB"
    ) %>%
    left_join(
      TFs[[gene]] %>%
        select(geneB, Rank_aggr_coexpr) %>%
        rename_with(~ paste0(., "_global"), c("Rank_aggr_coexpr")),
      by = "geneB"
    )
  
  # add literature curated TFs
  lit_TFs_filtered <- lit_TFs %>%
    filter(TF_Symbol == gene) %>%
    group_by(Target_Symbol) %>%
    summarise(
      Databases = paste(unique(Databases), collapse = ", "),
      TF_Species = paste(unique(TF_Species), collapse = ", "),
      Target_Species = paste(unique(Target_Species), collapse = ", "),
      Cell_Type = paste(unique(Cell_Type), collapse = ", "),
      Experiment_Type = paste(unique(Experiment_Type), collapse = ", "),
      PubMed_ID = paste(unique(PubMed_ID), collapse = ", ")
    ) %>%
    rename(geneB = Target_Symbol)
  
  # Merge literature TF data
  gene_merged <- gene_merged %>%
    left_join(lit_TFs_filtered, by = "geneB")
  
  # Apply delta threshold if provided
  if (!is.null(delta_threshold)) {
    gene_merged <- gene_merged %>% filter(delta_abs >= delta_threshold)
  }
  
  if (!is.null(mic_neuroinflammation)) {
    mic_neuroinflammation$geneB = mic_neuroinflammation$Gene
    mic_neuroinflammation$mic_neuroinflammation = mic_neuroinflammation$Function
    gene_merged <- gene_merged %>%
      left_join(mic_neuroinflammation %>% select(c(mic_neuroinflammation, geneB)), by = "geneB")
  }
  
  if (!is.null(SEA_AD_micGenes)) {
    SEA_AD_micGenes$SEA_AD_micGenes = SEA_AD_micGenes$GeneSet
    gene_merged <- gene_merged %>% 
      left_join(SEA_AD_micGenes %>% select(c(geneB, SEA_AD_micGenes)), by = "geneB")
  }
  
  if (!is.null(eRegulon) & gene %in% eRegulon$TF) {
    eRegulon$geneB = eRegulon$Gene 
    eRegulon$SCENIC_TF2G_rho = eRegulon$TF2G_rho
    eRegulon = eRegulon %>% filter(TF == gene) %>%
      select(c(SCENIC_TF2G_rho, geneB))
    
    gene_merged <- gene_merged %>%
      left_join(eRegulon, by = "geneB")
  }
  
  
  return(unique(gene_merged))
}





##### extract corr  ===================================================================================
### below is a complement to the above - it add cor values per study per condition 
filter_and_integrate_ranking <- function(delta_threshold = NULL, 
                                         ctl_threshold = NULL, 
                                         gene_name, 
                                         combined_tf_allTrg, 
                                         RUNX1_merged_all) {
  
  # Step 1: Apply filtering based on delta_threshold and ctl_threshold if specified
  RUNX1_merged_top <- RUNX1_merged_all
  
  if (!is.null(delta_threshold)) {
    RUNX1_merged_top <- RUNX1_merged_top %>% filter(delta_abs >= delta_threshold)
  }
  
  if (!is.null(ctl_threshold)) {
    RUNX1_merged_top <- RUNX1_merged_top %>% filter(freq_CTL >= ctl_threshold)
  }
  
  # Step 2: Filter the combined_tf_allTrg dataset for the specified gene
  gene_integrated_ranking_top_cor <- combined_tf_allTrg %>%
    filter(geneA == gene_name) %>%
    filter(geneB %in% RUNX1_merged_top$geneB) %>%
    left_join(
      RUNX1_merged_top %>%
        select(
          geneA, geneB, coExpr_rank_FZ, coExpr_rank_allRank, Rank_aggr_coexpr_global, 
          Databases, TF_Species, Target_Species, Cell_Type, Experiment_Type, PubMed_ID
        ),
      by = c("geneA", "geneB")
    )
  
  # Step 3: Clean the study column
  gene_integrated_ranking_top_cor$study <- gsub("(_[^_]+)$", "", gene_integrated_ranking_top_cor$study)
  
  return(gene_integrated_ranking_top_cor)
}



integrate_freq_ranking <- function(delta_threshold = NULL, 
                                         ctl_threshold = NULL, 
                                         gene_name, 
                                         combined_tf_allTrg, 
                                         RUNX1_merged_all) {
  
  # Step 1: Apply filtering based on delta_threshold and ctl_threshold if specified
  RUNX1_merged_top <- RUNX1_merged_all
  
  if (!is.null(delta_threshold)) {
    RUNX1_merged_top <- RUNX1_merged_top %>% filter(delta_abs >= delta_threshold)
  }
  
  if (!is.null(ctl_threshold)) {
    RUNX1_merged_top <- RUNX1_merged_top %>% filter(freq_CTL >= ctl_threshold)
  }
  
  # Step 2: Filter the combined_tf_allTrg dataset for the specified gene
  gene_integrated_ranking_top_cor <- combined_tf_allTrg %>%
    filter(geneA == gene_name) %>%
    filter(geneB %in% RUNX1_merged_top$geneB) %>%
    left_join(
      RUNX1_merged_top %>%
        select(
          geneA, geneB, freq_CTL, freq_AD
        ),
      by = c("geneA", "geneB")
    )
  
  # Step 3: Clean the study column
  gene_integrated_ranking_top_cor$study <- gsub("(_[^_]+)$", "", gene_integrated_ranking_top_cor$study)
  
  return(gene_integrated_ranking_top_cor)
}




calculate_differences <- function(data) {
  data %>%
    group_by(study, geneA, geneB) %>%
    summarise(
      delta_freq = unique(freq_CTL - freq_AD),
      delta_rank_norm_corr = diff(rank_norm_corr[order(condition)]),  # Ensure correct order (control, AD)
      .groups = "drop"
    )
}






























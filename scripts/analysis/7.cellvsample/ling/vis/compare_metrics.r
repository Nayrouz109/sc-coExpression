
library(data.table)



dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/ling-PEER4"
cells = c("Ast", "Inh", "Exc")

studies <- c("Lau-2020-ALL", "Lim-2022-ALL", 
  "ling-CTL", 
  "Nagy-2020-ALL",
  "Pineda-2021-ALL", "ROSMAP-Mathys-2019-ALL", "Velmeshev-2019-ALL"
)


read_and_tag <- function(study, source) {
  message("Reading ", study, " — ", source)
  file <- file.path(dir, sprintf("%s_PEER4_gene_correlations_with_background_%s.csv", study, source))
  fread(file) %>%
    filter(cell_type %in% cells) %>%
    mutate(source = source)
}


combined <- map(studies, function(study) {
  map_dfr(c("xCell", "xCellxSubject", "xSubject"), ~ read_and_tag(study, .x))
}) %>%
  set_names(studies)



combined1 <- lapply(combined, function(df) {
  setDT(df)
  
  # Compute background stats per source × cell_type
  bg_stats <- df[gene_list == "background_random1000",
                 .(bg_mean = mean(correlation, na.rm = TRUE),
                   bg_sd   = sd(correlation, na.rm = TRUE)),
                 by = .(source, cell_type)]
  
  # Merge stats back into main dataframe
  df <- merge(df, bg_stats, by = c("source", "cell_type"), all.x = TRUE)
  
  # Compute corrected correlations
  df[, corrected_cor_mean := correlation - bg_mean]
  df[, corrected_cor_z := (correlation - bg_mean) / bg_sd]
  
  # Rescale corrected correlations to [0, 1] within each source × cell_type
  #df[, corrected_cor_z_scaled :=
  #     if (all(is.na(corrected_cor_z))) NA_real_
  #   else (corrected_cor_z - min(corrected_cor_z, na.rm = TRUE)) /
  #     (max(corrected_cor_z, na.rm = TRUE) - min(corrected_cor_z, na.rm = TRUE)),
  #   by = .(source, cell_type)]
  
  # Optionally clean up
  df[, c("bg_mean", "bg_sd") := NULL]
  
  return(df)
})








###  compare studies  ==================================================================================== 
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)


### Step 1: Compute pairwise similarity metrics ==================================================================================== 
compare_studies <- function(df) {
  
  # Format wide study-edge matrix
  edge_matrix <- df %>%
    mutate(edge = paste(gene1, gene2, sep = "_")) %>%
    select(edge, study, corrected_cor_z) %>%
    pivot_wider(names_from = study, values_from = corrected_cor_z)
  
  studies <- colnames(edge_matrix)[-1]
  
  # Must have ≥2 studies
  if (length(studies) < 2) return(NULL)
  
  study_pairs <- combn(studies, 2, simplify = FALSE)
  
  map_dfr(study_pairs, function(p) {
    s1 <- p[1]
    s2 <- p[2]
    
    vals <- edge_matrix %>% select(all_of(c(s1, s2)))
    vals <- vals %>% filter(!is.na(.data[[s1]]) & !is.na(.data[[s2]]))
    
    pearson_cor <- suppressWarnings(cor(vals[[s1]], vals[[s2]], method = "pearson"))
    spearman_cor <- suppressWarnings(cor(vals[[s1]], vals[[s2]], method = "spearman"))
    
    pos_s1 <- which(vals[[s1]] > 0)
    pos_s2 <- which(vals[[s2]] > 0)
    jaccard_pos <- length(intersect(pos_s1, pos_s2)) / length(union(pos_s1, pos_s2))
    
    sign_preservation <- mean(sign(vals[[s1]]) == sign(vals[[s2]]))
    
    # --- Top-k overlap ---
    k = 200
    ord_s1 <- order(abs(vals[[s1]]), decreasing = TRUE)
    ord_s2 <- order(abs(vals[[s2]]), decreasing = TRUE)
    
    topk_s1 <- ord_s1[1:min(k, length(ord_s1))]
    topk_s2 <- ord_s2[1:min(k, length(ord_s2))]
    
    topk_overlap <- length(intersect(topk_s1, topk_s2)) / min(length(topk_s1), length(topk_s2))
    
    
    tibble(
      study1 = s1,
      study2 = s2,
      pearson_cor,
      spearman_cor,
      jaccard_pos,
      sign_preservation, 
      topk_overlap
    )
  })
}

### Step 3: Full workflow  ==================================================================================== 
run_one_combination <- function(df, source, cell_type, gene_list) {
  
  sub <- df %>%
    filter(
      source == !!source,
      cell_type == !!cell_type,
      gene_list == !!gene_list
    )
  
  res <- compare_studies(sub)
  if (is.null(res)) return(NULL)
  
  res %>%
    mutate(
      source = source,
      cell_type = cell_type,
      gene_list = gene_list,
      .before = study1
    )
}

run_all_comparisons <- function(df) {
  
  sources <- sort(unique(df$source))
  cell_types <- sort(unique(df$cell_type))
  gene_lists <- sort(unique(df$gene_list))
  
  combos <- expand_grid(
    source = sources,
    cell_type = cell_types,
    gene_list = gene_lists
  )
  
  map_dfr(seq_len(nrow(combos)), function(i) {
    run_one_combination(
      df,
      combos$source[i],
      combos$cell_type[i],
      combos$gene_list[i]
    )
  })
}


df <- dplyr::bind_rows(combined1, .id = "study")
results_all <- run_all_comparisons(df)














###  compare sources  ==================================================================================== 
compare_sources <- function(df) {
  
  # wide matrix: one column per source
  edge_matrix <- df %>%
    mutate(edge = paste(gene1, gene2, sep = "_")) %>%
    select(edge, source, corrected_cor_z) %>%
    pivot_wider(names_from = source, values_from = corrected_cor_z)
  
  sources <- colnames(edge_matrix)[-1]   # all source columns
  
  # Need ≥2 sources
  if (length(sources) < 2) return(NULL)
  
  source_pairs <- combn(sources, 2, simplify = FALSE)
  
  map_dfr(source_pairs, function(p) {
    s1 <- p[1]
    s2 <- p[2]
    
    vals <- edge_matrix %>% select(all_of(c(s1, s2)))
    vals <- vals %>% filter(!is.na(.data[[s1]]) & !is.na(.data[[s2]]))
    
    # correlation similarity
    pearson_cor  <- suppressWarnings(cor(vals[[s1]], vals[[s2]], method = "pearson"))
    spearman_cor <- suppressWarnings(cor(vals[[s1]], vals[[s2]], method = "spearman"))
    
    # overlap of positive edges
    pos_s1 <- which(vals[[s1]] > 0)
    pos_s2 <- which(vals[[s2]] > 0)
    jaccard_pos <- length(intersect(pos_s1, pos_s2)) /
      length(union(pos_s1, pos_s2))
    
    # sign preservation
    sign_preservation <- mean(sign(vals[[s1]]) == sign(vals[[s2]]))
    
    # --- Top-k strongest edges overlap ---
    k <- 200
    
    ord_s1 <- order(abs(vals[[s1]]), decreasing = TRUE)
    ord_s2 <- order(abs(vals[[s2]]), decreasing = TRUE)
    
    topk_s1 <- ord_s1[1:min(k, length(ord_s1))]
    topk_s2 <- ord_s2[1:min(k, length(ord_s2))]
    
    topk_overlap <- length(intersect(topk_s1, topk_s2)) /
      min(length(topk_s1), length(topk_s2))
    
    tibble(
      source1 = s1,
      source2 = s2,
      pearson_cor,
      spearman_cor,
      jaccard_pos,
      sign_preservation,
      topk_overlap
    )
  })
}

run_one_combination_sources <- function(df, study, cell_type, gene_list) {
  
  sub <- df %>%
    filter(
      study == !!study,
      cell_type == !!cell_type,
      gene_list == !!gene_list
    )
  
  res <- compare_sources(sub)
  if (is.null(res)) return(NULL)
  
  res %>%
    mutate(
      study = study,
      cell_type = cell_type,
      gene_list = gene_list,
      .before = source1
    )
}

run_all_source_comparisons <- function(df) {
  
  studies    <- sort(unique(df$study))
  cell_types <- sort(unique(df$cell_type))
  gene_lists <- sort(unique(df$gene_list))
  
  combos <- expand_grid(
    study = studies,
    cell_type = cell_types,
    gene_list = gene_lists
  )
  
  map_dfr(seq_len(nrow(combos)), function(i) {
    run_one_combination_sources(
      df,
      combos$study[i],
      combos$cell_type[i],
      combos$gene_list[i]
    )
  })
}



results_sources <- run_all_source_comparisons(df)








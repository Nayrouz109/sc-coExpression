


##### visualization ===================================================================================
#### correlation plots ===============================================================================
plot_delta_relationship <- function(data, title = "Gene Analysis") {
  ggplot(data, aes(x = delta_rank_norm_corr, y = delta_freq)) +
    geom_point(alpha = 0.3, size = 2) +  
    facet_wrap(~study, scales = "free") + 
    theme_minimal() +
    labs(
      title = title,
      x = "Delta Rank Norm Corr",
      y = "Delta Frequency"
    ) +
    theme(strip.text = element_text(size = 12, face = "bold"))
}

plot_delta_relationship <- function(data, 
                                    title = "Gene Analysis", 
                                    highlight_geneB = NULL) {
  ggplot(data, aes(x = delta_rank_norm_corr, y = delta_freq, color = (geneB == highlight_geneB))) +
    geom_point(alpha = 0.6, size = 2) +  
    facet_wrap(~study, scales = "free") + 
    theme_minimal() +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "firebrick"), guide = "none") + 
    labs(
      title = title,
      x = "Delta Rank Norm Corr",
      y = "Delta Frequency"
    ) +
    theme(strip.text = element_text(size = 12, face = "bold"))
}




plot_rank_norm_corr <- function(data, gene) {
  
  filtered_data <- data %>% filter(geneB == gene)
  
  # Create the plot
  p <- ggplot(filtered_data, aes(x = condition, y = rank_norm_corr)) +
    geom_violin(scale = "width", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_jitter(aes(color = study), width = 0.2, size = 2) +
    theme_classic() +
    labs(title = paste("Distribution of rank_norm_corr for", gene),
         x = "Condition",
         y = "Rank Normalized Correlation") +
    theme(legend.position = "right")
  
  return(p)
}



plot_rank_norm_corr <- function(data, gene) {
  
  
  filtered_data <- data %>% filter(geneB == gene)
  filtered_data$condition <- factor(filtered_data$condition, levels = c("control", "AD"))
  
  # Reshape data to wide format for paired t-test
  wide_data <- filtered_data %>%
    pivot_wider(names_from = condition, values_from = rank_norm_corr) %>%
    drop_na(control, AD)  # Ensure only pairs are considered
  
  # Perform paired t-test
  if (nrow(wide_data) > 1) {
    t_test_result <- t.test(wide_data$control, wide_data$AD, paired = TRUE)
    p_value <- signif(t_test_result$p.value, digits = 3)
  } else {
    p_value <- NA  # Not enough data for a paired test
  }
  
  # Create the plot
  p <- ggplot(filtered_data, aes(x = condition, y = rank_norm_corr, group = study)) +
    geom_point(size = 5, alpha = 0.7) +  # Increase point size
    geom_line(alpha = 0.5) +  # Connect points by study
    theme_classic(base_size = 15) +  # Increase text size
    labs(
      title = paste("Rank-Normalized Correlation for", gene),
      x = "Condition",
      y = "Rank Normalized Correlation",
      caption = ifelse(!is.na(p_value), paste("Paired t-test p-value:", p_value), "Not enough paired data for t-test")
    ) +
    theme(
      legend.position = "none",  # Remove legend
      axis.text = element_text(size = 14),  # Increase axis text size
      #axis.title = element_text(size = 22),  # Increase axis title size
      #plot.title = element_text(size = 26, face = "bold"),  # Increase title size
      plot.caption = element_text(size = 14, face = "italic")  # Increase caption size
    )
  
  return(p)
}
library(patchwork)










#### TF dot plots ===============================================================================
plot_tf_summary0 <- function(data, col1, col2) {
  tf_summary_df <- data %>% 
    mutate(gene_color = case_when(
      #!!sym(col1) & !!sym(col2) ~ "purple",  # If both conditions are true, color purple
      !!sym(col1) & !!sym(col2) ~ "firebrick",  # If both conditions are true, color purple
      !!sym(col1) ~ "firebrick",  # If only col1 is true, color red
      #!!sym(col2) ~ "blue",  # If only col2 is true, color blue
      !!sym(col2) ~ "firebrick",  # If only col2 is true, color blue
      TRUE ~ "black"))  # Default color is black
  
  summary_df <- head(tf_summary_df %>% 
                       select(geneA, meanTop200Ctl, gene_color) %>% 
                       arrange(desc(meanTop200Ctl)), 30)
  
  ggplot(tf_summary_df, aes(y = meanTop200Ctl, x = reorder(geneA, meanTop200Ctl))) +
    geom_point(aes(color = gene_color), size = 1) +  # Color points
    geom_text_repel(data = summary_df,
                    aes(x = geneA, y = meanTop200Ctl, label = geneA, fontface = "italic", color = gene_color), 
                    max.overlaps = 20,
                    force = 0.5,
                    nudge_x = -300,
                    hjust = 1,
                    direction = "y",
                    size = 5,
                    segment.size = 0.15,
                    segment.color = "black") +
    expand_limits(x = nrow(summary_df) + 50) +  # Prevent point cut off
    scale_color_identity() +  # Ensure colors are applied correctly
    ylab("mean Top200 in CTL") + 
    xlab("gene") + 
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
}



#### OLD OLD OLD ---------- OLD OLD OLD --------------- OLD OLD OLD -------------------------- #########
plot_tf_summary <- function(data, int_clmn = NULL) {
  if (!is.null(int_clmn) && int_clmn %in% colnames(data)) {
    tf_summary_df <- data %>% 
      mutate(gene_color = ifelse(!!sym(int_clmn) == TRUE, "firebrick", "black")) 
  } else {
    tf_summary_df <- data %>% mutate(gene_color = "black")
  }
  
  summary_df <- head(tf_summary_df %>% 
                       select(geneA, meanTop200Ctl, gene_color) %>% 
                       arrange(desc(meanTop200Ctl)), 30)
  
  ggplot(tf_summary_df, aes(y = meanTop200Ctl, x = reorder(geneA, meanTop200Ctl))) +
    geom_point(aes(color = gene_color), size = 1) +  # Color points
    geom_text_repel(data = summary_df,
                    aes(x = geneA, y = meanTop200Ctl, label = geneA, fontface = "italic", color = gene_color), 
                    max.overlaps = 20,
                    force = 0.5,
                    #nudge_x = -500,
                    nudge_x = -300,
                    hjust = 1,
                    direction = "y",
                    size = 5,
                    segment.size = 0.15,
                    segment.color = "black") +
    expand_limits(x = nrow(summary_df) + 50) +  # Prevent point cut off
    scale_color_identity() +  # Ensure colors are applied correctly
    ylab("mean Top200 in CTL") + 
    xlab("gene") + 
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
}

#### new below ---------- new below --------------- new below -------------------------- #########
plot_tf_summary <- function(data, null_data = NULL, int_clmn = NULL, cell_type = NULL) {
  
  # Assign colors based on int_clmn
  if (!is.null(int_clmn) && int_clmn %in% colnames(data)) {
    tf_summary_df <- data %>% 
      mutate(gene_color = ifelse(!!sym(int_clmn) == TRUE, "firebrick", "black")) 
  } else {
    tf_summary_df <- data %>% mutate(gene_color = "black")
  }
  
  # Select top 30 for labels
  summary_df <- head(tf_summary_df %>% 
                       select(geneA, meanTop200Ctl, gene_color) %>% 
                       arrange(desc(meanTop200Ctl)), 30)
  
  # Determine the best null overlap (e.g., maximum)
  #best_null_overlap <- if (!is.null(null_data)) {
  #  max(null_data$avg_overlap, na.rm = TRUE)
  #} else {
  #  NA
  #}
  # Add 2SD for the null data  
  null_data <- null_data %>%
    dplyr::mutate(
      null_mean = mean(avg_overlap, na.rm = TRUE),
      null_sd   = sd(avg_overlap, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      null_lo = null_mean - 2 * null_sd,
      null_hi = null_mean + 2 * null_sd
  )

  # Base plot
  p <- ggplot(tf_summary_df, aes(y = meanTop200Ctl, x = reorder(geneA, meanTop200Ctl))) +
    geom_point(aes(color = gene_color), size = 1) +
    geom_text_repel(data = summary_df,
                    aes(x = geneA, y = meanTop200Ctl, label = geneA, fontface = "italic", color = gene_color), 
                    max.overlaps = 20,
                    force = 0.5,
                    nudge_x = -300,
                    hjust = 1,
                    direction = "y",
                    size = 5,
                    segment.size = 0.15,
                    segment.color = "black") +
    expand_limits(x = nrow(summary_df) + 50) +
    scale_color_identity() +
    ylab("Mean Top200") + 
    xlab("Transcription Factor (TF)") + 
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
  
  # Add red horizontal line if null data is provided
  if (!is.null(null_data)) {
    #p <- p + geom_hline(yintercept = best_null_overlap, color = "red", linewidth = 1)
    # --- NEW: shaded band for null ±2SD ---
    p <- p +  annotate("rect",
                 xmin = -Inf, xmax = Inf,
                 ymin = null_summary$null_lo, ymax = null_summary$null_hi,
                 fill = "red", alpha = 0.1) +
    # --- NEW: null mean line ---
    geom_hline(yintercept = null_summary$null_mean, color = "red") +
    # --- NEW: null +2SD line ---
    geom_hline(yintercept = null_summary$null_hi, color = "red", linetype = "dashed")
    
  }
  
  # Add cell type annotation if provided
  if (!is.null(cell_type)) {
    p <- p + annotate("text", 
                      x = max(as.numeric(factor(tf_summary_df$geneA))), 
                      y = min(tf_summary_df$meanTop200Ctl, na.rm = TRUE) - 10, 
                      label = cell_type, 
                      hjust = 1, vjust = 0, 
                      size = 8, fontface = "bold")
  }
  
  return(p)
}







### TF, ranks heatmap ==================================================================================
library(tidyr)
library(stringr)
library(purrr)
library(pheatmap)
library(tibble)
library(viridis)

plot_tf_heatmap <- function(data_list, top_n = 30, viridis_option = "D") {
  
  tf_elements <- grep("rpr_tf", names(data_list), value = TRUE)
  
  tf_df <- purrr::map_dfr(tf_elements, function(name) {
    df <- data_list[[name]]
    cell_type <- stringr::str_extract(name, "^[^_]+")
    df %>%
      dplyr::select(geneA, meanTop200Ctl, brain_TF) %>%
      dplyr::mutate(cell_type = cell_type)
  })
  
  heatmap_data <- tf_df %>%
    tidyr::pivot_wider(names_from = cell_type, values_from = meanTop200Ctl)
  
  heatmap_matrix <- heatmap_data %>%
    tibble::column_to_rownames("geneA") %>%
    as.matrix()
  
  ranked_matrix <- apply(heatmap_matrix[, -1], 2, rank, ties.method = "average")
  
  top_genes_list <- apply(ranked_matrix, 2, function(x) {
    names(sort(x, decreasing = TRUE))[1:top_n]
  })
  
  top_genes_union <- unique(as.vector(top_genes_list))
  
  subset_matrix <- ranked_matrix[rownames(ranked_matrix) %in% top_genes_union, ]
  
  brain_tf_annot <- tf_df %>%
    dplyr::distinct(geneA, brain_TF) %>%
    tibble::column_to_rownames("geneA")
  
  row_annot <- data.frame(brain_TF = brain_tf_annot[rownames(subset_matrix), , drop = FALSE])
  row_annot$brain_TF <- as.character(row_annot$brain_TF)
  
  ann_colors <- list(brain_TF = c("TRUE" = "firebrick", "FALSE" = "black"))
  
  viridis_colors <- viridis::viridis(100, option = viridis_option)
  
  pheatmap::pheatmap(
    subset_matrix, 
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    treeheight_row = 0,
    treeheight_col = 0, 
    annotation_row = row_annot,
    annotation_colors = ann_colors,
    color = viridis_colors,
    fontsize_row = 6,
    fontsize_col = 10,
    fontface_col = "bold",         # ✅ Makes cell type labels bold
    main = paste("Top", top_n, "TFs across cell types"), 
    annotation_names_row = FALSE,  # ✅ Hides the "brain_TF" label bar
    border_color = NA
  )
}





library(dplyr)
library(ggplot2)
library(forcats)
library(purrr)

plot_lollipop_tf_frequency <- function(data_list, top_n = 30, dot_size = 4, line_thickness = 1) {
  # Filter only "rpr_tf" elements
  tf_dfs <- data_list[grep("rpr_tf", names(data_list))]
  
  # Extract cell type from the name
  extract_cell_type <- function(name) {
    strsplit(name, "_")[[1]][1]
  }
  
  # Combine into one long ranked dataframe
  tf_summary_df <- imap_dfr(tf_dfs, function(df, name) {
    cell_type <- extract_cell_type(name)
    df %>%
      mutate(
        geneA = as.character(geneA),
        meanTop200_AD = as.numeric(meanTop200_AD)
      ) %>%
      arrange(desc(meanTop200_AD)) %>%
      mutate(rank = row_number(), cell_type = cell_type) %>%
      select(geneA, rank, cell_type)
  })
  
  # Filter to top N per cell type
  top_tf_df <- tf_summary_df %>% filter(rank <= top_n)
  
  # Count number of cell types each TF appears in
  tf_count <- top_tf_df %>%
    group_by(geneA) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))
  
  # Lollipop plot (TFs on x-axis)
  ggplot(tf_count, aes(x = fct_reorder(geneA, count), y = count)) +
    geom_segment(aes(xend = geneA, y = 0, yend = count), color = "black", size = line_thickness, linetype = "dotted") +
    geom_point(color = "darkblue", size = dot_size) +
    labs(
      title = paste0("TF occurrence among top ", top_n, " most reproducible across cell types"),
      x = "Transcription Factor (TF)",
      y = paste0("Number of cell types") #with top ", top_n, " reproducibility")
    ) +
    #theme_bw(base_size = 14) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) + 
    theme(#axis.text = element_text(size = 20),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
}



plot_lollipop_tf_frequency_colored <- function(data_list, top_n = 30, dot_size = 4, line_thickness = 1) {
  # Filter only "rpr_tf" elements
  tf_dfs <- data_list[grep("rpr_tf", names(data_list))]
  
  # Extract cell type from the name
  extract_cell_type <- function(name) {
    strsplit(name, "_")[[1]][1]
  }
  
  # Combine into one long ranked dataframe
  tf_summary_df <- imap_dfr(tf_dfs, function(df, name) {
    cell_type <- extract_cell_type(name)
    df %>%
      mutate(
        geneA = as.character(geneA),
        meanTop200_AD = as.numeric(meanTop200_AD)
      ) %>%
      arrange(desc(meanTop200_AD)) %>%
      mutate(rank = row_number(), cell_type = cell_type) %>%
      select(geneA, rank, cell_type)
  })
  
  # Filter to top N per cell type
  top_tf_df <- tf_summary_df %>% filter(rank <= top_n)
  
  # Count number of cell types each TF appears in
  tf_count <- top_tf_df %>%
    group_by(geneA) %>%
    summarise(count = n(),
              cell_types = paste(cell_type, collapse = ","),  # track where it occurs
              .groups = "drop") %>%
    arrange(desc(count))
  
  # --- Custom ordering ---
  tf_count <- tf_count %>%
    mutate(ordering_var = ifelse(count == 1, cell_types, paste0("zzz_", count))) %>%
    arrange(ordering_var, desc(count)) %>%
    mutate(geneA = factor(geneA, levels = unique(geneA)))
  
  # --- Build band info for count=1 genes ---
  band_df <- tf_count %>%
    filter(count == 1) %>%
    group_by(cell_types) %>%
    summarise(
      xmin = min(as.numeric(geneA)) - 0.5,
      xmax = max(as.numeric(geneA)) + 0.5,
      ypos = max(count),   # height of band = max count in this group (here always 1)
      xpos = mean(as.numeric(geneA)),
      .groups = "drop"
    )
  
  # --- Plot ---
  ggplot(tf_count, aes(x = geneA, y = count)) +
    # Bands that stop at lollipop height
    geom_rect(
      data = band_df,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ypos + 0.3, fill = cell_types),
      inherit.aes = FALSE, alpha = 0.3
    ) +
    # Lollipop lines
    geom_segment(aes(xend = geneA, y = 0, yend = count),
                 color = "black", size = line_thickness, linetype = "dotted") +
    # Dots for count=1 (colored by band, solid alpha)
    geom_point(
      data = tf_count %>% filter(count == 1),
      aes(color = cell_types),
      size = dot_size, alpha = 1
    ) +
    # Dots for count>1 (always black)
    geom_point(
      data = tf_count %>% filter(count > 1),
      color = "black",
      size = dot_size
    ) +
    # Labels above bands
    geom_text(
      data = band_df,
      aes(x = xpos, y = ypos + 0.6, label = cell_types),
      inherit.aes = FALSE, fontface = "bold", size = 5
    ) +
    labs(
      title = paste0("TF occurrence among top ", top_n, " most reproducible across cell types"),
      x = "Transcription Factor (TF)",
      y = "Number of cell types"
    ) +
    scale_fill_brewer(palette = "Set3") +
    scale_color_brewer(palette = "Set3") +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 16), 
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20),
      legend.position = "none",
      plot.margin = margin(c(20, 10, 10, 10))
    )
}


##### upset plot =========================================
library(UpSetR)
#gene_sets <- split(top_tf_df$geneA, top_tf_df$cell_type)
#text.scale = c(
  #1.5, # intersection size (top numbers above bars)
  #1.5, # intersection labels (x-axis below matrix)
  #1.5, # set size (left numbers beside set bars)
  #1.5, # set names (left labels)
  #1.2, # main title
  #1.2  # axis title
#)

#upset(
  #fromList(gene_sets),
  #sets = names(gene_sets),
  #order.by = "freq",
  #mainbar.y.label = "Shared Genes",
  #point.size = 5,
  #line.size = 0.8,
  #text.scale = c(1.8, 1.8, 1.8, 1.8, 1.5, 1.5),
  #keep.order = TRUE,
  #show.numbers = "yes",
  #number.angles = 0,
  ##set_size.show = FALSE,
  ##sets.bar.color = NA,
  #att.pos = "none"
#)












summarize_top_tfs <- function(data_list, top_n = 30) {
  # Filter only "rpr_tf" elements
  tf_dfs <- data_list[grep("rpr_tf", names(data_list))]
  
  # Extract cell type from the name
  extract_cell_type <- function(name) {
    strsplit(name, "_")[[1]][1]
  }
  
  # Combine into one long ranked dataframe
  tf_summary_df <- imap_dfr(tf_dfs, function(df, name) {
    cell_type <- extract_cell_type(name)
    
    df %>%
      mutate(
        geneA = as.character(geneA),  # Ensure geneA is a character
        meanTop200_AD = as.numeric(meanTop200_AD)  # Ensure numeric
      ) %>%
      arrange(desc(meanTop200_AD)) %>%
      mutate(
        rank = row_number(),
        cell_type = cell_type
      ) %>%
      select(geneA, rank, cell_type)
  })
  
  # Filter to top N per cell type
  top_tf_df <- tf_summary_df %>%
    filter(rank <= top_n)
  
  # Count number of cell types each TF appears in
  tf_count <- top_tf_df %>%
    group_by(geneA) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))
  
  # Plot barplot of TF frequency across cell types
  # Plot with TFs on x-axis (no flipping)
  ggplot(tf_count, aes(x = fct_reorder(geneA, count), y = count)) +
    geom_col(fill = "steelblue") +
    labs(
      title = paste0("TF occurrence among top ", top_n, " most reproducible across cell types"),
      x = "TF (geneA)",
      y = paste0("Number of cell types with top ", top_n, " reproducibility")
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1)
    )
  
}








plot_tf_rank_vs_expr <- function(tf_rpr, tf_xpr, top_genes_union = NULL, label_genes = NULL) {
  
  # Ensure tf_rpr is a data.frame
  tf_rpr <- as.data.frame(tf_rpr)
  
  # Rank and scale within cell type
  tf_rpr <- tf_rpr %>%
    group_by(cell_type) %>%
    mutate(
      mean_ranked_rpr = rank(meanTop200Ctl, ties.method = "average"),
      mean_ranked_rpr_scaled = (mean_ranked_rpr - min(mean_ranked_rpr)) / 
        (max(mean_ranked_rpr) - min(mean_ranked_rpr))
    ) %>%
    ungroup()
  
  # Average expression per TF per cell type
  tf_xpr_avg_by_celltype <- tf_xpr[, .(
    mean_expr_log_cpm = mean(avg_expr_log_cpm, na.rm = TRUE)
  ), by = .(geneA, cell_type)]
  
  # Merge
  merged_tf <- merge(
    tf_xpr_avg_by_celltype,
    tf_rpr,
    by = c("geneA", "cell_type"),
    all = FALSE
  )
  
  merged_tf$brain_TF <- as.logical(merged_tf$brain_TF)
  merged_tf$TF <- paste0(merged_tf$geneA, "_", merged_tf$cell_type)
  
  # Defaults
  if (is.null(top_genes_union)) top_genes_union <- character(0)
  if (is.null(label_genes)) label_genes <- character(0)

  # Plot
  p <- ggplot(as.data.frame(merged_tf), aes(x = mean_expr_log_cpm, y = mean_ranked_rpr_scaled)) +
    geom_point(aes(
      #color = ifelse(geneA %in% top_genes_union, "firebrick", "black"),
      #size = ifelse(geneA %in% top_genes_union, 2, 1),
      #alpha = ifelse(geneA %in% top_genes_union, 0.8, 0.2),
      #color = ifelse(geneA %in% label_genes, "firebrick", "black"),
      fill = ifelse(geneA %in% label_genes, "firebrick", "white"), 
      size = ifelse(geneA %in% label_genes, 4, 1.5),
      alpha = ifelse(geneA %in% label_genes, 0.9, 0.4)
    ),
    color = "black",  # Border color
    shape = 21, stroke = 0.5) +  # Filled circle
    #shape = 1, stroke = 1) +
    
    # Label selected genes
    geom_text_repel(
      data = subset(merged_tf, geneA %in% label_genes),
      aes(label = TF), size = 5, box.padding = 0.5, max.overlaps = 30
    ) +
    
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      legend.position = "none"
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_identity() +
    scale_alpha_identity()
  
  return(p)
}



### TF, ribo, null ==========================================================================================
library(ggplot2)
library(dplyr)

plot_TF_ribo_null <- function(data_list, selected_elements, null_name, pos_ctrl_name) {
  # Combine selected elements into one dataframe
  combined_df <- bind_rows(
    lapply(selected_elements, function(name) {
      df <- data_list[[name]]
      df$source <- name  # Track source element
      
      # Categorize as 'null', 'positive_control', or 'actual'
      df$data_type <- ifelse(name == null_name, "Null",
                             ifelse(name == pos_ctrl_name, "Positive Control", "Actual Data"))
      return(df)
    })
  )
  
  ggplot(combined_df, aes(x = source, y = meanTop200Ctl, fill = source)) +
    #geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
    geom_jitter(aes(color = data_type), width = 0.35, alpha = 0.6, size = 2.5) +  
    labs(x = "List Elements", y = "meanTop200Ctl", title = "Comparison of meanTop200Ctl") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x labels
    scale_fill_brewer(palette = "Set2") +  # Boxplot colors
    scale_color_manual(values = c("Null" = "gray50", "Positive Control" = "orange", "Actual Data" = "black")) +   # Custom jitter colors
    theme_classic() +
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20),
          legend.position = "none")
}

plot_TF_ribo_null_Pct <- function(data_list, selected_elements, null_name, pos_ctrl_name) {
  # Combine selected elements into one dataframe
    combined_df <- bind_rows(
    lapply(selected_elements, function(name) {
      df <- data_list[[name]]
      df$source <- name  # Track source element
      
      # Categorize as 'null', 'positive_control', or 'actual'
      df$data_type <- ifelse(name == null_name, "Null",
                             ifelse(name == pos_ctrl_name, "Ribo", "TR"))
      return(df)
    })
  )
  ribo_n <- combined_df %>% filter(data_type == "Ribo") %>% nrow()

  # Create a new column with conditional logic
  combined_df <- combined_df %>%
    mutate(meanTop200Pct = case_when(
      data_type %in% c("Null", "TR") ~ (meanTop200Ctl / 200) * 100,
      data_type == "Ribo" ~ (meanTop200Ctl / ribo_n) * 100,
      TRUE ~ NA_real_
  ))
  
  ggplot(combined_df, aes(x = data_type, y = meanTop200Pct, fill = data_type)) +
    geom_jitter(aes(color = data_type), width = 0.35, alpha = 0.6, size = 2.5) +  
    labs( y = "% Mean of Top Partners", 
          #title = "Comparison of Normalized meanTop200Ctl"
          ) +
    theme_classic() +
    theme(axis.text.x = element_text(#angle = 0, 
                                      #hjust = 1
                                      ),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        legend.position = "none") +
    scale_fill_brewer(palette = "Set2") +
    scale_color_manual(values = c("Null" = "gray50", "Ribo" = "orange", "TR" = "black"))
}


# Example usage:
#selected_elements <- c("Mic_tf_reprod_tbl", "Ast_tf_reprod_tbl", "Exc_tf_reprod_tbl")
#plot_TF_ribo_null(data_list, selected_elements, null_name = "Ast_tf_reprod_tbl", pos_ctrl_name = "Exc_tf_reprod_tbl")


#### TF expression variance ===============================================================================

# 1. variance of TF expression across datasets per cell type 

plot_tf_rank_vs_expr_variance <- function(tf_rpr, tf_xpr, top_genes_union = NULL, label_genes = NULL) {
  library(data.table) 
  
  # Ensure tf_rpr is a data.frame
  tf_rpr <- as.data.frame(tf_rpr)
  
  # Rank and scale reproducibility scores within each cell type
  tf_rpr <- tf_rpr %>%
    group_by(cell_type) %>%
    mutate(
      mean_ranked_rpr = rank(meanTop200Ctl, ties.method = "average"),
      mean_ranked_rpr_scaled = (mean_ranked_rpr - min(mean_ranked_rpr)) / 
        (max(mean_ranked_rpr) - min(mean_ranked_rpr))
    ) %>%
    ungroup()
  
  # Calculate expression variance per gene (TF) per cell type
  tf_xpr_variance_by_celltype <- tf_xpr[, .(
    expr_variance = var(avg_expr_log_cpm, na.rm = TRUE)
  ), by = .(geneA, cell_type)]
  
  # Merge reproducibility with expression variance
  merged_tf <- merge(
    tf_xpr_variance_by_celltype,
    tf_rpr,
    by = c("geneA", "cell_type"),
    all = FALSE
  )
  
  merged_tf$brain_TF <- as.logical(merged_tf$brain_TF)
  merged_tf$TF <- paste0(merged_tf$geneA, "_", merged_tf$cell_type)
  
  # Defaults for genes to highlight
  if (is.null(top_genes_union)) top_genes_union <- character(0)
  if (is.null(label_genes)) label_genes <- character(0)
  
  # Plot expression variance (x) vs reproducibility rank (y)
  p <- ggplot(as.data.frame(merged_tf), aes(x = expr_variance, y = mean_ranked_rpr_scaled)) +
    geom_point(aes(
      fill = ifelse(geneA %in% label_genes, "firebrick", "white"),
      size = ifelse(geneA %in% label_genes, 4, 1.5),
      alpha = ifelse(geneA %in% label_genes, 0.9, 0.4)
    ),
    color = "black", shape = 21, stroke = 0.5) +
    
    geom_text_repel(
      data = subset(merged_tf, geneA %in% label_genes),
      aes(label = TF), size = 5, box.padding = 0.5, max.overlaps = 30
    ) + 
    
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "none"
    ) +
    labs(
      x = "Expression Variance (log CPM)",
      y = "Scaled Reproducibility Rank"
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_identity() +
    scale_alpha_identity()
  
  return(p)
}


# 2. variance of TF expression across cell types
plot_tf_rpr_vs_cross_celltype_expr_var <- function(tf_rpr, tf_xpr, label_genes = NULL) {
  
  tf_rpr <- as.data.frame(tf_rpr)
  
  # Rank and scale reproducibility scores within each cell type
  tf_rpr <- tf_rpr %>%
    group_by(cell_type) %>%
    mutate(
      mean_ranked_rpr = rank(meanTop200Ctl, ties.method = "average"),
      mean_ranked_rpr_scaled = (mean_ranked_rpr - min(mean_ranked_rpr)) /
        (max(mean_ranked_rpr) - min(mean_ranked_rpr))
    ) %>%
    ungroup()
  
  # Calculate expression mean per gene per cell type (if not already averaged)
  dt <- as.data.table(tf_xpr)
  expr_avg_by_gene_celltype <- dt[, .(
    mean_expr = mean(avg_expr_log_cpm, na.rm = TRUE)
  ), by = .(geneA, cell_type)]
  
  # Calculate expression variance across cell types per gene
  expr_var_across_celltypes <- expr_avg_by_gene_celltype[, .(
    expr_variance_across_celltypes = var(mean_expr, na.rm = TRUE)
  ), by = geneA]
  
  # Merge variance back with reproducibility per gene per cell type
  merged <- merge(tf_rpr, expr_var_across_celltypes, by = "geneA", all.x = TRUE)
  
  merged$TF <- paste0(merged$geneA, "_", merged$cell_type)
  
  # Plot reproducibility vs cross-cell-type expression variance
  p <- ggplot(merged, aes(x = expr_variance_across_celltypes, y = mean_ranked_rpr_scaled)) +
    geom_point(aes(
      fill = ifelse(geneA %in% label_genes, "firebrick", "white"),
      size = ifelse(geneA %in% label_genes, 4, 1.5),
      alpha = ifelse(geneA %in% label_genes, 0.9, 0.4)
    ),
    color = "black", shape = 21, stroke = 0.5) +
    
    geom_text_repel(
      data = subset(merged, geneA %in% label_genes),
      aes(label = TF), size = 5, box.padding = 0.5, max.overlaps = 30
    ) +
    
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.position = "none"
    ) +
    labs(
      x = "Expression Variance Across Cell Types (log CPM)",
      y = "Scaled Reproducibility Rank (within Cell Type)"
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_identity() +
    scale_alpha_identity()
  
  return(p)
}

#### TF reproducibaility plots ===============================================================================

plot_tf_reprod <- function(data, int_genes = NULL) {
  ggplot(data, aes(x = meanTop200_AD, y = meanTop200Ctl, color = ifelse(geneA %in% int_genes, "highlight", "default"))) +
    geom_point(alpha = 0.4, size = 3) +  
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +  # Solid red y = x line
    scale_color_manual(values = c("highlight" = "blue", "default" = "black")) +
    theme_classic() +
    ylab("mean Top200 in CTL") + 
    xlab("mean Top200 in AD") + 
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20),
          legend.position = "none")
}

# plot specific columns
plot_tf_reprod <- function(data, int_genes = NULL, x_col, y_col) {
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], 
                   color = ifelse(geneA %in% int_genes, "highlight", "default"))) +
    geom_point(alpha = 0.4, size = 3) +  
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +  
    scale_color_manual(values = c("highlight" = "blue", "default" = "black")) +
    theme_classic() +
    ylab(paste("Mean Top200 in", y_col)) + 
    xlab(paste("Mean Top200 in", x_col)) + 
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20),
          legend.position = "none")
}

# adds labels 
plot_tf_reprod2 <- function(data, int_genes = NULL, label_genes = NULL) {
  library(ggrepel)
  ggplot(data, aes(x = meanTop200_AD, y = meanTop200Ctl, color = ifelse(geneA %in% int_genes, "highlight", "default"))) +
    geom_point(alpha = 0.4, size = 3) +  
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +  # Solid red y = x line
    scale_color_manual(values = c("highlight" = "blue", "default" = "black")) +
    geom_text_repel(data = subset(data, geneA %in% label_genes),
                    aes(label = geneA), size = 6, box.padding = 0.5, max.overlaps = Inf) +
    theme_classic() +
    ylab("mean Top200 in CTL") + 
    xlab("mean Top200 in AD") + 
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 20),
          legend.position = "none")
}


#### unique targets per TF vs TF reproducability ==============================================================================================================
## 1. what kind of cor is between TF rpr score and number of Top200 union per cell type 
plot_tf_rpr_union <- function(df, interesting_genes = NULL) {
  
  # Flag interesting genes
  df <- df %>%
    mutate(
      is_interesting = if (!is.null(interesting_genes)) geneA %in% interesting_genes else FALSE,
      color = ifelse(is_interesting, "firebrick", "black"),
      size = ifelse(is_interesting, 4, 2.5),
      alpha = ifelse(is_interesting, 0.9, 0.4),
      label = ifelse(is_interesting, paste0(geneA, "_", cell_type), NA)
    )
  
  # Plot
  p <- ggplot(df, aes(x = n_geneB, y = tf_rpr_score)) +
    # Black points first
    geom_point(
      data = subset(df, !is_interesting),
      aes(color = color, size = size, alpha = alpha)
    ) +
    # Red points on top
    geom_point(
      data = subset(df, is_interesting),
      aes(color = color, size = size, alpha = alpha)
    ) +
    ggrepel::geom_label_repel(
      data = subset(df, is_interesting),
      aes(label = label),
      size = 5,
      color = "black",
      fill = "white",
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.color = "black"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text = element_text(size = 15, face = "bold"),
      axis.title = element_text(size = 16),
      legend.position = "none"
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_size_identity() +
    scale_alpha_identity()
  
  return(p)
}




## 2. cell type TF-gene frequency 
plot_tf_pairfreq <- function(df, tf_of_interest, n_dataset_limit = 6, alpha_col = "avg_rank_norm_corr") {
  
  # Filter the data
  tf_data <- df %>%
    filter(geneA == tf_of_interest, n_datasets >= n_dataset_limit)
  
  # Create plot
  p <- ggplot(tf_data, aes(x = cell_type, y = pair_freq)) +
    geom_violin(fill = "grey80", color = NA, alpha = 0.8) +
    geom_jitter(
      aes(
        size = binding_score,
        alpha = .data[[alpha_col]]  # dynamic alpha column
      ),
      width = 0.1,
      shape = 21,
      fill = "firebrick",
      color = "black",
      stroke = 0.3
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    xlab("Cell Type") +
    ylab("Pair Frequency") +
    scale_size_continuous(name = "Binding Score") +
    scale_alpha_continuous(range = c(0.3, 0.9), name = alpha_col) +
    annotate(
      "text",
      x = 0.5,  # fractional position: leftmost
      y = max(tf_data$pair_freq, na.rm = TRUE) * 1.125,  # slightly above top
      label = tf_of_interest,
      hjust = 0,  # left-aligned
      vjust = 1,  # top-aligned
      size = 6,
      fontface = "bold"
    )
  
  return(p)
}






plot_tf_pairfreq <- function(df, tf_of_interest, n_dataset_limit = 6) {
  
  # Filter the data
  tf_data <- df %>%
    filter(geneA == tf_of_interest, n_datasets >= n_dataset_limit)
  
  # Create plot
  p <- ggplot(tf_data, aes(x = cell_type, y = pair_freq)) +
    geom_violin(fill = "grey80", color = NA, alpha = 0.8) +
    geom_jitter(
      aes(size = binding_score, alpha = avg_rank_norm_corr),
      width = 0.1,
      shape = 21,
      fill = "firebrick",
      color = "black",
      stroke = 0.3
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    xlab("Cell Type") +
    ylab("Pair Frequency") +
    scale_size_continuous(name = "Binding Score") +
    scale_alpha_continuous(range = c(0.3, 0.9), name = "Avg Rank Norm Corr") +
    annotate(
      "text",
      x = 0.5,  # fractional position: leftmost
      y = max(tf_data$pair_freq, na.rm = TRUE) * 1.125,  # slightly above top
      label = tf_of_interest,
      hjust = 0,  # left-aligned
      vjust = 1,  # top-aligned
      size = 6,
      fontface = "bold"
    )
  
  return(p)
}


## 3. unique targets
plot_tf_target_vs_rpr <- function(tf_target_dt, geneA_interest, 
                                  dataset_threshold = NULL, pair_freq_threshold = NULL) {
  
  # apply conditional filtering
  filtered_dt <- tf_target_dt
  if (!is.null(dataset_threshold)) {
    filtered_dt <- filtered_dt %>% filter(n_datasets >= dataset_threshold)
  }
  if (!is.null(pair_freq_threshold)) {
    filtered_dt <- filtered_dt %>% filter(pair_freq == pair_freq_threshold)
  }
  
  # summarize
  unique_targets_summary <- filtered_dt %>%
    group_by(geneA, cell_type) %>%
    summarize(
      n_uniqueTargets = n_distinct(geneB),
      tf_rpr_score = first(tf_rpr_score),
      .groups = "drop"
    )
  
  # prepare styling
  unique_targets_summary <- unique_targets_summary %>%
    mutate(
      highlight = ifelse(geneA == geneA_interest, "firebrick", "grey80"),
      pt_size   = ifelse(geneA == geneA_interest, 4, 2),
      pt_alpha  = ifelse(geneA == geneA_interest, 0.9, 0.4)
    )
  
  # plot
  ggplot(unique_targets_summary, aes(x = n_uniqueTargets, y =  tf_rpr_score)) +
    geom_point(
      data = subset(unique_targets_summary, geneA != geneA_interest),
      aes(size = pt_size, alpha = pt_alpha),
      fill = "grey80", shape = 21, stroke = 0.5
    ) +
    geom_point(
      data = subset(unique_targets_summary, geneA == geneA_interest),
      aes(size = pt_size, alpha = pt_alpha),
      fill = "firebrick", shape = 21, stroke = 0.5
    ) +
    geom_label_repel(
      data = subset(unique_targets_summary, geneA == geneA_interest),
      aes(label = paste(geneA, cell_type)),
      size = 5, box.padding = 0.5, max.overlaps = 30,
      color = "black", fill = "white"
    ) +
    labs(
      title = paste0("TF reproducibility vs. unique targets: ", geneA_interest),
      x = "Number of unique coexpressed targets (n_uniqueTargets)",
      y =  "TF reproducibility score (tf_rpr_score)"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text   = element_text(size = 15, face = "bold"),
      axis.title  = element_text(size = 16),
      legend.position = "none"
    ) +
    scale_size_identity() +
    scale_alpha_identity()
}




#### recent =================================================================================
plot_pairwise_cor <- function(cor_list, geneA_interest, geneB_interest) {
  library(dplyr)
  library(ggplot2)
  library(purrr)
  
  # Step 1: extract correlations for geneA-geneB pair across cell types
  df <- map_dfr(names(cor_list), function(ct) {
    cor_df <- cor_list[[ct]]
    
    # filter for both possible orientations of the pair
    cor_pair <- cor_df %>%
      filter(
        (geneA == geneA_interest & geneB == geneB_interest) |
          (geneA == geneB_interest & geneB == geneA_interest)
      ) %>%
      mutate(CellType = ct)
    
    return(cor_pair)
  })
  
  # Step 2: average correlation per cell type
  avg_df <- df %>%
    group_by(CellType) %>%
    summarise(avg_rank_norm_corr = mean(avg_rank_norm_corr, na.rm = TRUE), .groups = "drop")
  
  # Step 3: plot
  p <- ggplot(avg_df, aes(x = CellType, y = avg_rank_norm_corr, group = 1)) +
    geom_line(color = "firebrick", size = 1) +
    geom_point(color = "firebrick", size = 3) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text   = element_text(size = 15, face = "bold"),
      axis.title  = element_text(size = 16)
    ) +
    labs(
      title = paste0("Avg correlation of ", geneA_interest, " & ", geneB_interest),
      x = "Cell type",
      y = "Average correlation"
    )
  
  return(list(data = avg_df, plot = p))
}







#### 1 TF  ===============================================================================
#ggplot(RUNX1_merged_all, aes(x = freq_AD, y = freq_CTL)) +
#  geom_point(color = "black", alpha = 0.4, size = 3) +  
#  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +  # Solid red y = x line
#  theme_classic() +
#  ylab ("frequency of occurance in CTL") + 
#  xlab ("frequency of occurance in AD") + 
#  theme(axis.text = element_text(size = 20), 
#        axis.title = element_text(size = 20))
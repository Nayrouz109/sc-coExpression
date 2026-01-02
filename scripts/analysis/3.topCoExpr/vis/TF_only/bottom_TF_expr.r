
#### this is related to bottomw TFs expression visalization 
plot_tf_expression_by_celltype <- function(tf_xpr, genes_of_interest) {
  # Ensure data is a data.table
  tf_xpr <- as.data.table(tf_xpr)
  
  # Filter for multiple genes
  gene_data <- tf_xpr[geneA %in% genes_of_interest]
  
  # Plot
  p <- ggplot(gene_data, aes(
    x = cell_type,
    y = avg_expr_log_cpm,
    group = study_name,
    color = study_name
  )) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    #facet_wrap(~ geneA, scales = "free_y") +   # ← Each TF gets its own panel
    facet_wrap(~ geneA, scales = "fixed") +   # ← Each TF gets its own panel
    labs(
      title = "TF Expression Across Cell Types",
      x = "Cell Type",
      y = "Average Expression (log CPM)",
      color = "Study"
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  return(p)
}


plot_tf_expression_by_celltype(tf_xpr, c("ZNF282", "ZNF333", "ZNF317", "ZNF121"))
plot_tf_expression_by_celltype(tf_xpr, c("PRDM4", "PRDM10", "GMEB1", "GMEB2", "SETDB1", "HMGXB4", "MBD1"))
plot_tf_expression_by_celltype(tf_xpr, c("USF1", "USF3", "E4F1", "REL"))


### cell type specificity 
plot_tf_expression_by_celltype(tf_xpr, c("ESRRA", "NR1H3", "NR2C2", "NR6A1"))
plot_tf_expression_by_celltype(tf_xpr, c("OLIG1", "OLIG2", "SOX10", "POU6F1", "FOXP4"))
plot_tf_expression_by_celltype(tf_xpr, c("ESRRA", "NR1H3", "NR2C2", "NR6A1", "RUNX1", "STAT2"))
##########################################################################################################
##########################################################################################################
##########################################################################################################

plot_tf_expression_combined <- function(tf_xpr, genes_of_interest) {
  # Ensure data is a data.table
  tf_xpr <- as.data.table(tf_xpr)
  
  # Filter for the genes of interest
  gene_data <- tf_xpr[geneA %in% genes_of_interest]
  
  # Plot
  p <- ggplot(gene_data, aes(
    x = cell_type,
    y = avg_expr_log_cpm,
    group = interaction(geneA, study_name), # separate lines per study but same color per gene
    color = geneA
  )) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = "TF Expression Across Cell Types and Studies",
      x = "Cell Type",
      y = "Average Expression (log CPM)",
      color = "TF"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  return(p)
}
plot_tf_expression_combined(tf_xpr, c("ZNF282", "PRDM4", "PRDM10", "REL"))

















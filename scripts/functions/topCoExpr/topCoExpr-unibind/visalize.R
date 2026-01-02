
library(data.table)
library(ggplot2)


## 1.0 TF/gene per cell type ===================================================================================
plot_tf_expression_by_celltype <- function(tf_xpr, gene_of_interest) {
  # Ensure data is a data.table
  tf_xpr <- as.data.table(tf_xpr)
  
  # Filter for the gene of interest
  gene_data <- tf_xpr[geneA == gene_of_interest]
  
  # Plot
  p <- ggplot(gene_data, aes(x = cell_type, y = avg_expr_log_cpm, group = study_name, color = study_name)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = paste("Expression of", gene_of_interest, "across Cell Types"),
      x = "Cell Type",
      y = "Average Expression (log CPM)",
      color = "Study"
    ) +
    theme_minimal(base_size = 14) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20)
    )
  
  return(p)
}

## 1.1 TF-gene pair per cell type ===================================================================================

plot_two_genes_expression <- function(tf_xpr, genes_of_interest) {
  # Ensure tf_xpr is a data.table
  tf_xpr <- as.data.table(tf_xpr)
  
  # Ensure input is exactly 2 genes
  if (length(genes_of_interest) != 2) {
    stop("Please provide exactly two genes.")
  }

  # Assign colors
  gene_colors <- setNames(c("firebrick", "gray40"), genes_of_interest)
  #gene_colors <- setNames(c("darkorange2", "gray"), genes_of_interest)
  #gene_colors <- setNames(c("orangered", "chartreuse4"), genes_of_interest)
  
  # Filter data for the two genes
  gene_data <- tf_xpr[gene %in% genes_of_interest]
  
  # Plot
  p <- ggplot(gene_data, aes(x = cell_type, y = avg_expr_log_cpm,
                             group = interaction(gene, study_name),
                             color = gene
                             )) +
    geom_line(#aes(linetype = study_name), 
              size = 0.5) +
    geom_point(size = 2) +
    scale_color_manual(values = gene_colors) + 
    labs(
      title = paste("Expression of", paste(genes_of_interest, collapse = " & "), "across Cell Types"),
      x = "Cell Type",
      y = "Average Expression (log CPM)",
      color = "Gene",
      linetype = "Study"
    ) +
    theme_minimal(base_size = 14) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
  
  return(p)
}

## 2 TF-gene binding vs coexpression ===================================================================================
plot_binding_vs_coexpr <- function(df, geneA_query, geneB_query, filter_tf = TRUE) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  
  if (filter_tf) {
    df <- df %>% filter(geneA == geneA_query)
  }
  
  df_long <- df %>%
    pivot_longer(cols = Ast:Opc, names_to = "CellType", values_to = "CoexprReproducibility") %>%
    mutate(highlight = geneA == geneA_query & geneB == geneB_query)
  
  df_bg <- df_long %>% filter(!highlight)
  df_fg <- df_long %>% filter(highlight)
  
  ggplot() +
    # Background points
    geom_point(
      data = df_bg,
      aes(x = CoexprReproducibility, y = binding_score),
      shape = 21, fill = "lightgrey", color = "black", size = 3, alpha = 0.4, stroke = 0.2
    ) +
    # Highlighted points
    geom_point(
      data = df_fg,
      aes(x = CoexprReproducibility, y = binding_score),
      shape = 21, fill = "firebrick", color = "black", size = 5, alpha = 0.9, stroke = 0.2
    ) +
    # Add readable labels with repel
    geom_text_repel(
      data = df_fg,
      aes(x = CoexprReproducibility, y = binding_score, label = CellType),
      #color = "firebrick", 
      box.padding = 0.8, 
      size = 4, 
      #fontface = "bold", 
      #max.overlaps = Inf, 
      max.overlaps = 30
    ) +
    labs(
      title = paste("Binding Score vs Co-expression for", geneA_query, "-", geneB_query),
      x = "Co-expression Reproducibility",
      y = "Binding Score"
    ) +
    theme_minimal(base_size = 14)
}

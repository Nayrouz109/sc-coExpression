



### this is to visualize TF reproducability scores across different cell types 

source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/visalize.R")
### 0. needed data ===================================================================================================
dir <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell"
files <- list.files(dir, pattern = ".*rpr.*.rds", full.names = TRUE)

data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))

### TF top 30 per cell type ===================================================================================================
### ribo vs null ===================================================================================================
cell_types <- unique(sub("_t200.*", "", names(data_list)))
ribo_null <- list()

for (ct in cell_types) {
  selected_elements <- grep(paste0("^", ct, "_t200_rpr_"), names(data_list), value = TRUE)
  null_name <- paste0(ct, "_t200_rpr_tn")
  pos_ctrl_name <- paste0(ct, "_t200_rpr_rb")
  
  ribo_null[[ct]][["selected_elements"]] = selected_elements
  ribo_null[[ct]][["null_name"]] = null_name
  ribo_null[[ct]][["pos_ctrl_name"]] = pos_ctrl_name 
}

library(cowplot)
# Astrocytes
combined_ast <- ggdraw() + 
  draw_plot(plot_tf_summary(data_list$Ast_t200_rpr_tf, null_data = data_list$Ast_t200_rpr_tn, int_clmn = "brain_TF", cell_type = "Astrocytes")) + 
  draw_plot(plot_TF_ribo_null_Pct(data_list, selected_elements = ribo_null$Ast$selected_elements, null_name = ribo_null$Ast$null_name, pos_ctrl_name = ribo_null$Ast$pos_ctrl_name), x = 0.13, y = 0.45, width = 0.4, height = 0.5)

# Excitatory neurons
combined_exc <- ggdraw() + 
  draw_plot(plot_tf_summary(data_list$Exc_t200_rpr_tf, null_data = data_list$Exc_t200_rpr_tn, int_clmn = "brain_TF", cell_type = "Excitatory neurons")) + 
  draw_plot(plot_TF_ribo_null_Pct(data_list, selected_elements = ribo_null$Exc$selected_elements, null_name = ribo_null$Exc$null_name, pos_ctrl_name = ribo_null$Exc$pos_ctrl_name), x = 0.13, y = 0.45, width = 0.4, height = 0.5)

# Inhibitory neurons
combined_inh <- ggdraw() + 
  draw_plot(plot_tf_summary(data_list$Inh_t200_rpr_tf, null_data = data_list$Inh_t200_rpr_tn, int_clmn = "brain_TF", cell_type = "Inhibitory neurons")) + 
  draw_plot(plot_TF_ribo_null_Pct(data_list, selected_elements = ribo_null$Inh$selected_elements, null_name = ribo_null$Inh$null_name, pos_ctrl_name = ribo_null$Inh$pos_ctrl_name), x = 0.13, y = 0.45, width = 0.4, height = 0.5)

# Microglia
combined_mic <- ggdraw() + 
  draw_plot(plot_tf_summary(data_list$Mic_t200_rpr_tf, null_data = data_list$Mic_t200_rpr_tn, int_clmn = "brain_TF", cell_type = "Microglia")) + 
  draw_plot(plot_TF_ribo_null_Pct(data_list, selected_elements = ribo_null$Mic$selected_elements, null_name = ribo_null$Mic$null_name, pos_ctrl_name = ribo_null$Mic$pos_ctrl_name), x = 0.13, y = 0.45, width = 0.4, height = 0.5)

# Oligodendrocytes
combined_oli <- ggdraw() + 
  draw_plot(plot_tf_summary(data_list$Oli_t200_rpr_tf, null_data = data_list$Oli_t200_rpr_tn, int_clmn = "brain_TF", cell_type = "Oligodendrocytes")) + 
  draw_plot(plot_TF_ribo_null_Pct(data_list, selected_elements = ribo_null$Oli$selected_elements, null_name = ribo_null$Oli$null_name, pos_ctrl_name = ribo_null$Oli$pos_ctrl_name), x = 0.13, y = 0.45, width = 0.4, height = 0.5)

# OPCs
combined_opc <- ggdraw() + 
  draw_plot(plot_tf_summary(data_list$Opc_t200_rpr_tf, null_data = data_list$Opc_t200_rpr_tn, int_clmn = "brain_TF", cell_type = "OPCs")) + 
  draw_plot(plot_TF_ribo_null_Pct(data_list, selected_elements = ribo_null$Opc$selected_elements, null_name = ribo_null$Opc$null_name, pos_ctrl_name = ribo_null$Opc$pos_ctrl_name), x = 0.13, y = 0.45, width = 0.4, height = 0.5)


combined_ast
combined_exc
combined_inh
combined_mic
combined_oli
combined_opc

#728 by 900



### TF rpr heatmap + lollipop ===================================================================================================
viridis_options <- tibble::tibble(
  option = c("A", "B", "C", "D", "E", "F", "G", "H"),
  name = c("magma", "inferno", "plasma", "viridis", "cividis", "rocket", "mako", "turbo")
)

plot_tf_heatmap(data_list, top_n = 30, viridis_option = "E") #290 by 900 

plot_lollipop_tf_frequency(data_list, top_n = 30, dot_size = 4, line_thickness = 0.5)
plot_lollipop_tf_frequency_colored(data_list, top_n = 30, dot_size = 4, line_thickness = 0.5) #1700 by 600
plot_lollipop_tf_frequency_colored(data_list, top_n = 50, dot_size = 4, line_thickness = 0.5) #1900 by 600
plot_lollipop_tf_frequency_colored(data_list, top_n = 100, dot_size = 4, line_thickness = 0.5) #3000 by 800

#generated internally from the plot_tf_heatmap
saveRDS(tf_df, "/home/nelazzabi/rewiring/data/coExpr/0all-human/TF_rpr_inegrated_cellTypes.rds")
saveRDS(top_genes_union, "/home/nelazzabi/rewiring/data/coExpr/0all-human/TF_rpr_tp30_genes_union.rds")




### TF expr vs rpr score ===================================================================================================

tf_rpr = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell/intgr/TF_rpr_inegrated_cellTypes.rds") %>% setDT()
tf_xpr = fread("/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/TF_exp_summary_logCPM.csv") %>% setDT()
top_genes_union = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/TF_rpr_tp30_genes_union.rds")

tf_celltype_dt <- data.table(
  TF = c(
    "TBR1", "NEUROD1", "NEUROD2", "NPAS4", "MEF2C", "NEUROG2",         # Excitatory
    "ASCL1", "LHX6", "DLX1", "DLX2", "GATA2", "GATA3", "NR4A1",        # Inhibitory
    "SPI1", "SALL1", "IRF8", "MEF2C", "KLF4"                           # Microglia
  ),
  CellType = c(
    rep("Excitatory Neuron", 6),
    rep("Inhibitory Neuron", 7),
    "Microglia", "Microglia", "Microglia", "Microglia", "Microglia"
  ))

#tf_celltype_dt$TF[(tf_celltype_dt$TF %in% tf_rpr$geneA) ]
# "MEF2C" "NR4A1" "MEF2C"https://rstudio.pavlab.msl.ubc.ca/graphics/plot_zoom_png?width=694&height=878

plot_tf_rank_vs_expr(tf_rpr , 
                     tf_xpr , 
                     top_genes_union, label_genes = c("STAT4"))


plot_tf_rank_vs_expr(tf_rpr , 
                     tf_xpr , 
                     top_genes_union, label_genes = NULL)

plot_tf_rank_vs_expr(tf_rpr , 
                     tf_xpr , 
                     top_genes_union, label_genes = c("NPAS3", "RORA", "NCOR2"))


### TF expr per cell type ===================================================================================================
library(ggplot2)
library(data.table)

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

plot_tf_expression_by_celltype(tf_xpr, "RUNX1")
plot_tf_expression_by_celltype(tf_xpr, "STAT2")
plot_tf_expression_by_celltype(tf_xpr, "NFIB")

plot_tf_expression_by_celltype(tf_gene_xpr, "ARHGAP26")



# TF expression variance  ==============================================================================================
# 1. variance of TF expression across datasets per cell type   
plot_tf_rank_vs_expr_variance (tf_rpr, tf_xpr, top_genes_union = NULL, label_genes = NULL) 
plot_tf_rank_vs_expr_variance (tf_rpr, tf_xpr, top_genes_union = NULL, label_genes = c("RUNX1", "STAT2")) 
plot_tf_rank_vs_expr_variance (tf_rpr, tf_xpr, top_genes_union = NULL, label_genes = c("RUNX1")) 
plot_tf_rank_vs_expr_variance (tf_rpr, tf_xpr, top_genes_union = NULL, label_genes = c("STAT2")) 

# 2. variance of TF expression across cell types
plot_tf_rpr_vs_cross_celltype_expr_var(tf_rpr, tf_xpr, label_genes = c("RUNX1", "STAT2"))
plot_tf_rpr_vs_cross_celltype_expr_var(tf_rpr, tf_xpr, label_genes = c("ZFP36L1"))








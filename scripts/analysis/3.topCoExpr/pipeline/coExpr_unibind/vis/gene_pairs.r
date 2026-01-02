

source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/topCoExpr-unibind/visalize.R")


tf_gene_xpr = fread("/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/allGenes_exp_summary_logCPM.csv") %>% setDT()
tf_top_with_binding = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human-unibind/TF_top_with_binding_ranked_scaled.rds")


#### 1. general trend ====================================================================
df_long <- tf_top_with_binding %>%
  pivot_longer(cols = Ast:Opc, names_to = "CellType", values_to = "CoexprReproducibility") 

plot(df_long$binding_score, df_long$CoexprReproducibility) 


### 2. plot gene pairs  ==================================================================
plot_two_genes_expression(tf_gene_xpr, c("RUNX1", "ARHGAP26"))
plot_tf_expression_by_celltype(tf_gene_xpr %>% mutate(geneA = gene), "RUNX1")

plot_binding_vs_coexpr(tf_top_with_binding, "RUNX1", "ARHGAP26", filter_tf = TRUE)


plot_two_genes_expression(tf_gene_xpr, c("RUNX1", "GRB2"))
plot_binding_vs_coexpr(tf_top_with_binding, "RUNX1", "GRB2", filter_tf = TRUE)


plot_two_genes_expression(tf_gene_xpr, c("RUNX1", "B3GNT5"))
plot_binding_vs_coexpr(tf_top_with_binding, "RUNX1", "B3GNT5", filter_tf = TRUE)


plot_two_genes_expression(tf_gene_xpr, c("RUNX1", "TBC1D14"))
plot_binding_vs_coexpr(tf_top_with_binding, "RUNX1", "TBC1D14", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("RUNX1", "ARHGAP12"))
plot_binding_vs_coexpr(tf_top_with_binding, "RUNX1", "ARHGAP12", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("RUNX1", "RNF130"))
plot_binding_vs_coexpr(tf_top_with_binding, "RUNX1", "RNF130", filter_tf = TRUE)


plot_two_genes_expression(tf_gene_xpr, c("SP1", "RGS12"))
plot_binding_vs_coexpr(tf_top_with_binding, "SP1", "RGS12", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("JUN", "ENC1"))
plot_binding_vs_coexpr(tf_top_with_binding, "JUN", "ENC1", filter_tf = TRUE)

### TFs  ========================================================================
View(tf_top_with_binding %>% filter(geneA == "RORA"))
View(tf_top_with_binding %>% filter(geneA == "FOXN3"))
View(tf_top_with_binding %>% filter(geneA == "FOXP2"))
View(tf_top_with_binding %>% filter(geneA == "NFIB"))

### NFIB  ========================================================================

plot_two_genes_expression(tf_gene_xpr, c("NFIB", "SCD5")) #* 
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "SCD5", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("NFIB", "HDAC4")) #* 
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "HDAC4", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("NFIB", "HOMER2")) 
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "HOMER2", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("NFIB", "PDE7A")) 
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "PDE7A", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("NFIB", "NFIA"))
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "NFIA", filter_tf = TRUE)


plot_two_genes_expression(tf_gene_xpr, c("NFIB", "DUSP16"))
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "DUSP16", filter_tf = TRUE)

### DGE + little messy 
plot_two_genes_expression(tf_gene_xpr, c("NFIB", "GPR37"))
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "GPR37", filter_tf = TRUE)

### less than 10 
plot_two_genes_expression(tf_gene_xpr, c("NFIB", "ATP2B4"))
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "ATP2B4", filter_tf = TRUE)

plot_two_genes_expression(tf_gene_xpr, c("NFIB", "USP31"))
plot_binding_vs_coexpr(tf_top_with_binding, "NFIB", "USP31", filter_tf = TRUE)








### STAT2 ========================================================================
## there are patterns of differential gene expression    = oli vs mic 
## there are patterns of no differential gene expression = oli vs exc 

## TF expression remains the comparable 
## target changes in expression 
plot_two_genes_expression(tf_gene_xpr, c("STAT2", "FAM171A1"))
plot_binding_vs_coexpr(tf_top_with_binding, "STAT2", "FAM171A1", filter_tf = TRUE)


# no binding evidence 
plot_two_genes_expression(tf_gene_xpr, c("STAT2", "RASGRF2"))
plot_binding_vs_coexpr(tf_top_with_binding, "STAT2", "RASGRF2", filter_tf = TRUE)


# DE of target in all other cell types 
plot_two_genes_expression(tf_gene_xpr, c("STAT2", "LRP2"))
plot_binding_vs_coexpr(tf_top_with_binding, "STAT2", "LRP2", filter_tf = TRUE)







####### TF-gene target overlap between cells ==================================================================================
dir <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell"
files <- list.files(dir, pattern = ".*rpr.*.rds", full.names = TRUE)
data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))

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
top_tf_df <- tf_summary_df %>% filter(rank <= 30)
top30TF_union = unique(top_tf_df$geneA) 


####### TF target overlap between cells ==================================================================================

library(dplyr)
int_gene = top30TF_union
int_gene = "NFIA" 

top_per_celltype <- function(cell_type) {
  tf_top_with_binding %>% filter(geneA %in% int_gene) %>%
    filter(.data[[cell_type]] > 0) %>%         # Keep only pairs active in this cell type
    arrange(desc(.data[[cell_type]])) %>%      # Sort by that cell type's value
    slice_head(n = 200) %>%
    select(gene_pair) %>%
    mutate(cell_type = cell_type)
}

cell_types <- c("Ast", "Exc", "Inh", "Mic", "Oli", "Opc")
top200_list <- lapply(cell_types, top_per_celltype)
top200_df <- bind_rows(top200_list)

library(tidyr)
overlap_matrix <- top200_df %>%
  pivot_wider(names_from = cell_type, values_from = cell_type,
              values_fn = length, values_fill = 0) %>%
  mutate(across(-gene_pair, ~ ifelse(. > 0, 1, 0)))  # Ensure binary



### heatmap 
library(pheatmap)
library(dplyr)
library(tibble)

# Prepare matrix
mat <- overlap_matrix %>%
  column_to_rownames("gene_pair") %>%
  as.matrix()

# Define custom color palette
custom_colors <- c("0" = "black", "1" = "blue")

# Convert to numeric to factor for color mapping
pheatmap(mat,
         color = c("black", "blue"),
         breaks = c(-0.01, 0.5, 1.01),  # ensures correct binning of 0 and 1
         show_rownames = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = int_gene, 
         legend_breaks = c(0,1),
         legend_labels = c("Absent", "Present"), 
         fontsize_col = 20)


### unique counts ----******* -------- ****** ---------


# Step 1: Convert matrix to tibble and keep rownames
mat_df <- mat %>%
  as_tibble(rownames = "gene_pair")

# Step 2: Identify rows where gene pair appears in only one cell type (sum == 1)
unique_rows <- mat_df %>%
  rowwise() %>%
  mutate(num_cell_types = sum(c_across(Ast:Opc) == 1)) %>%
  ungroup() %>%
  filter(num_cell_types == 1) %>%
  select(-num_cell_types)

unique_counts <- unique_rows %>%
  pivot_longer(cols = Ast:Opc, names_to = "cell_type", values_to = "value") %>%
  filter(value == 1) %>%
  group_by(cell_type) %>%
  summarise(n_unique_gene_pairs = n()) %>%
  arrange(desc(n_unique_gene_pairs))

ggplot(unique_counts, aes(x = reorder(cell_type, -n_unique_gene_pairs), y = n_unique_gene_pairs)) +
  geom_col(fill = "black") +
  geom_text(aes(label = n_unique_gene_pairs), vjust = -0.3, fontface = "bold") +
  labs(
    title = "Number of Gene Pairs Unique to Each Cell Type",
    x = "Cell Type",
    y = "Unique Gene Pairs"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 15), 
        axis.text.y = element_text(face = "bold" , size = 15) )



#### unique targets per TF vs TF reproducability ==============================================================================================================
tf_target_dt = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/intg_TF_200Trg_avgRnkCoexpr_binding_lit_TFrpr.rds")

unique_targets = tf_target_dt %>% filter(pair_freq == 1) #%>% group_by(geneA) %>% slice_head(n = 200)
unique_targets_summary <- unique_targets %>%
  group_by(geneA, cell_type) %>%
  summarize(
    n_uniqueTargets = n_distinct(geneB),
    tf_rpr_score = unique(tf_rpr_score),
    .groups = "drop"
  )







plot_tf_target_vs_rpr <- function(summary_df, geneA_interest) {
  summary_df$highlight <- ifelse(summary_df$geneA == geneA_interest, "firebrick", "white")
  summary_df$pt_size <- ifelse(summary_df$geneA == geneA_interest, 4, 1.5)
  summary_df$pt_alpha <- ifelse(summary_df$geneA == geneA_interest, 0.9, 0.4)
  
  ggplot(summary_df, aes(x = tf_rpr_score, y = n_uniqueTargets)) +
    geom_point(
      aes(
        fill = highlight,
        size = pt_size,
        alpha = pt_alpha
      ),
      color = "black", shape = 21, stroke = 0.5
    ) +
    geom_text_repel(
      data = subset(summary_df, geneA == geneA_interest),
      aes(label = paste(geneA, cell_type)),
      size = 5, box.padding = 0.5, max.overlaps = 30
    ) +
    labs(
      title = paste0("TF reproducibility vs. unique targets: ", geneA_interest),
      x = "TF reproducibility score (tf_rpr_score)",
      y = "Number of unique coexpressed targets (n_uniqueTargets)"
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
}

plot_tf_target_vs_rpr(unique_targets_summary, geneA_interest = "RUNX1")









### shared counts 
cell_types <- unique(top200_df$cell_type)
shared_counts <- combn(cell_types, 2, function(x) {
  a <- top200_df %>% filter(cell_type == x[1]) %>% pull(gene_pair)
  b <- top200_df %>% filter(cell_type == x[2]) %>% pull(gene_pair)
  tibble(cell_type_1 = x[1], cell_type_2 = x[2], shared = length(intersect(a, b)))
}, simplify = FALSE) %>%
  bind_rows()

shared_counts <- shared_counts %>%
  mutate(pair = paste(cell_type_1, cell_type_2, sep = " vs "))

ggplot(shared_counts, aes(x = reorder(pair, -shared), y = shared)) +
  geom_col(fill = "blue") +
  geom_text(aes(label = shared), vjust = -0.3, size = 5, fontface = "bold") +
  labs(
    title = "Number of Shared Top 200 Gene Pairs Between Cell Types",
    x = "Cell Type Pair",
    y = "Number of Shared Gene Pairs"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 15), 
        axis.text.y = element_text(face = "bold" , size = 15) )



























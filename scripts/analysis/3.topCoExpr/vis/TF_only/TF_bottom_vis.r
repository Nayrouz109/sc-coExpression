source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/visalize.R")
### 0. needed data ===================================================================================================
dir <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell"
files <- list.files(dir, pattern = ".*rpr.*.rds", full.names = TRUE)
data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))
data_list <- data_list[grep("rpr_tf", names(data_list))]


tf_long <- imap_dfr(data_list, ~
                      mutate(.x,
                             CellType = str_remove(.y, "_t200_rpr_tf"))) %>% select(-c("meanTop200_AD")) %>%  
  group_by(CellType) %>%
  mutate(
    Rank = rank(meanTop200Ctl, ties.method = "min")) %>%
  ungroup() 
  
get_low_tf_summary <- function(tf_long, quantile_cutoff = 0.75) {
  tf_long %>%
    #group_by(CellType) %>%
    #mutate(is_low = Rank <= quantile(Rank, quantile_cutoff)) %>%
    #ungroup() %>%
    group_by(geneA) %>%
    summarise(
      #n_celltypes_low   = sum(is_low),
      #n_celltypes_total = n_distinct(CellType),
      #frac_low          = n_celltypes_low / n_celltypes_total,
      avg_Rank          = mean(Rank, na.rm = TRUE)
    ) #%>%
  #mutate(
  #Category = case_when(
  #frac_low > 0.5 ~ "Globally lowly ranked TF",
  #frac_low > 0 & frac_low < 0.5    ~ "Cell-type-specific lowly ranked TF",
  #TRUE ~ "Not lowly ranked"
  #)
  #) %>%
  #arrange(desc(frac_low)) #%>%
  #filter(Category != "Not lowly ranked")
}


low_tf_summary <- get_low_tf_summary(tf_long, quantile_cutoff = 0.75)





  






##### plot lowly ranked TFs =======================================================================
### 700 by 1200 
library(ggplot2)
library(ggrepel)
library(dplyr)

### bottem avg ranked TFs 
bottom_df <- low_tf_summary %>%
  arrange(avg_Rank) %>%
  head(500) %>%
  mutate(group = "bottom") %>% 
  filter(avg_Rank < 400)
# i copied and pasted the above 500 into chatGPT and asked to categorize them 
# Combined TF examples

tf_examples1 <- c(
  # Zinc-finger proteins
  "ZNF384", "ZNF557", "ZNF514", "ZNF282", "ZNF333",
  "ZNF778", "ZNF383", "ZNF561", "ZNF143", "ZNF317",
  "ZNF805", "ZNF623", "ZNF583", "ZSCAN21", "ZSCAN12",
  "ZBTB17", "ZBTB26", "ZNF140", "ZNF445", "ZNF250", 
  # Phosphorylation-regulated TFs
  "ATF1", "SP2", "E2F6", 
  # Epigenetic regulators
  "PIAS3", "MBD1", "SETDB1"
)

tf_examples2 <- c(
  # Epigenetic regulators
   "SETDB2", "MBD2", "PIAS2",
  "PIAS4", "ZBTB4",
  
  # Phosphorylation-regulated TFs
  "CREB1", "CREM", "ATF2", "ATF6", "ATF6B",
  "E2F4", "SP1", "ELK1", "ELK4",
  "STAT1", "STAT2", "STAT5B", "FOXP4", "SRF", "MEF2D",
  
  # Hormone/signal-responsive TFs
  "NR1H2", "NR1H3", "NR2C1", "NR2C2", "ARNT"
)

right_df = low_tf_summary %>% filter(low_tf_summary$geneA %in% tf_examples1) %>% mutate(group = "right")
left_df = low_tf_summary %>% filter(low_tf_summary$geneA %in% tf_examples2) %>% mutate(group = "left")

label_df <- bind_rows(right_df, left_df) %>% mutate(nudge_x_val = ifelse(group == "left", -500, 650))

# Single plot
p <- ggplot(low_tf_summary, aes(y = avg_Rank, x = reorder(geneA, avg_Rank))) +
  geom_point(size = 1) +
  geom_text_repel(
    data = label_df,
    aes(x = geneA, y = avg_Rank, label = geneA, fontface = "italic"),
    nudge_x = label_df$nudge_x_val,
    nudge_y = 0.5, 
    hjust = 1,
    direction = "y",
    force = 1,
    max.iter = 10000,
    max.overlaps = 90,
    size = 5,
    segment.size = 0.15,
    segment.color = "black"
  ) +
  scale_color_identity() +
  ylab("Mean Top200") + 
  xlab("Transcription Factor (TF)") + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20),
    plot.margin = margin(c(10, 10, 10, 10))
  )

p



#### 000. heatmap example of Low vs High cell vs High Global  -------------------------------------------------------
all_low = c(tf_examples1, tf_examples2)
high_global = (global_summary %>% filter(Category == "Global core"))$geneA # refer back to TF_TOP_vis2 
high_cell <- global_summary %>% filter(Category == "Cell-type specific", !grepl("^(ZN|ZK|ZS)", geneA)) %>% pull(geneA) # refer back to TF_TOP_vis2 

add_genes <- c("ATF2", "ATF6", "CREM", "NR2C1")
remove_genes <- c("AFF1", "AFF2", "AFF3", "AHR", "ARNT2", "ARID5B", "ARID1B")
high_cell <- high_cell %>% union(add_genes) %>% # add (without creating duplicates)
  setdiff(remove_genes)      # remove unwanted genes

# identify cell type of max rank for the "High in one cell type" genes
cell_specific_order <- tf_long %>%
  filter(geneA %in% high_cell) %>%
  group_by(geneA) %>%
  slice_max(order_by = Rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(CellType, desc(Rank)) %>%   # grouped by cell type, then highest rank
  pull(geneA)

# order for the other two blocks
low_order   <- sort(c(tf_examples1, "CARF", "GABPB1"))
global_order <- sort(high_global)
# final ordering
all_gene_order <- c(low_order, global_order, cell_specific_order)

tf_long2 <- tf_long %>%
  mutate(Category_group = case_when(
    geneA %in% low_order ~ "Low in all cell types",
    geneA %in% high_global ~ "High in all cell types",
    geneA %in% high_cell ~ "High in one cell type"
  )) %>%
  filter(geneA %in% all_gene_order) %>%
  mutate(
    Category_group = factor(Category_group,
                            levels = c("Low in all cell types",
                                       "High in all cell types",
                                       "High in one cell type")),
    geneA = factor(geneA, levels = all_gene_order)
  )

ggplot(tf_long2, aes(x = geneA, y = CellType, fill = Rank)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno") + #cividis
  facet_grid(. ~ Category_group, scales = "free_x", space = "free_x") +
  scale_y_discrete(position = "right") +   # <- move CellType labels to the right
  labs(
    x = "Transcription Factor",
    y = "Cell Type",
    fill = "Norm. Reproducibility"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.background = element_rect(fill = "white"),
    #strip.text.x = element_text(size = 11, face = "bold"), 
    strip.text.x = element_blank() 
  )

# 1500 by 250 








### deleted plot ==============================================================
# Top 30 high frac_low TFs
top_df <- low_tf_summary %>%
  arrange(desc(frac_low)) %>%
  head(149) %>%
  mutate(group = "top")

# Bottom 30 low frac_low TFs
bottom_df <- low_tf_summary %>%
  arrange(frac_low) %>%
  head(200) %>%
  mutate(group = "bottom")

# Combine for labeling
label_df <- bind_rows(top_df, bottom_df) %>%
  mutate(nudge_x_val = ifelse(group == "top", -350, 350))


# Single plot
p <- ggplot(low_tf_summary, aes(y = avg_Rank, x = reorder(geneA, avg_Rank))) +
  geom_point(size = 1) +
  geom_text_repel(
    data = label_df,
    aes(x = geneA, y = avg_Rank, label = geneA, fontface = "italic"),
    nudge_x = label_df$nudge_x_val,
    hjust = 1,
    direction = "y",
    force = 1,
    max.iter = 10000,
    max.overlaps = 90,
    size = 5,
    segment.size = 0.15,
    segment.color = "black"
  ) +
  scale_color_identity() +
  ylab("Mean Top200") + 
  xlab("Transcription Factor (TF)") + 
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20),
    plot.margin = margin(c(10, 10, 10, 10))
  )

p



### xpr plots ==============================================================
tf_rpr = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell/intgr/TF_rpr_inegrated_cellTypes.rds") %>% setDT()
tf_xpr = fread("/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/TF_exp_summary_logCPM.csv") %>% setDT()
top_genes_union = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/TF_rpr_tp30_genes_union.rds")



plot_tf_rank_vs_expr <- function(tf_rpr, tf_xpr, label_genes = NULL) {
  
  # ensure df
  tf_rpr <- as.data.frame(tf_rpr)
  
  # rank + scale
  tf_rpr <- tf_rpr %>%
    group_by(cell_type) %>%
    mutate(mean_ranked_rpr = rank(meanTop200Ctl, ties.method = "average"),
           mean_ranked_rpr_scaled = (mean_ranked_rpr - min(mean_ranked_rpr)) / 
             (max(mean_ranked_rpr) - min(mean_ranked_rpr))) %>%
    ungroup()
  
  # avg expression per TF per cell type
  tf_xpr_avg <- tf_xpr[, .(
    mean_expr_log_cpm = mean(avg_expr_log_cpm, na.rm = TRUE)
  ), by = .(geneA, cell_type)]
  
  # merge
  merged_tf <- merge(tf_xpr_avg, tf_rpr, by = c("geneA", "cell_type"))
  merged_tf$TF <- paste0(merged_tf$geneA, "_", merged_tf$cell_type)
  
  # defaults
  if (is.null(label_genes)) label_genes <- character(0)
  
  # color palette for multiple label genes
  label_pal <- setNames(
    RColorBrewer::brewer.pal(max(3, length(label_genes)), "Dark2")[seq_along(label_genes)],
    label_genes
  )
  
  
  okabe_dark <- c(
    "#0072B2", "#D55E00", "#009E73", "green" , "blue", "orange" ,"black", "purple", 
    "#E69F00", "#CC79A7", "#56B4E9", "firebrick"
  )
  
  label_pal <- setNames(
    okabe_dark[seq_len(length(label_genes))],
    label_genes
  )
  
  
  merged_tf$label_color <- ifelse(
    merged_tf$geneA %in% label_genes,
    label_pal[merged_tf$geneA],
    "grey60"
  )
  
  # point size and alpha
  merged_tf$label_size <- ifelse(merged_tf$geneA %in% label_genes, 4.5, 2)
  merged_tf$label_alpha <- ifelse(merged_tf$geneA %in% label_genes, 1, 0.4)
  
  # unique shape per cell type
  merged_tf$cell_type <- factor(merged_tf$cell_type)
  
  p <- ggplot(merged_tf, aes(
    x = mean_expr_log_cpm,
    y = mean_ranked_rpr_scaled,
    shape = cell_type
  )) +
    geom_point(
      aes(
        color = geneA,
        size  = label_size,
        alpha = label_alpha
      )
    ) +
    scale_color_manual(
      values = c(
        "other" = "grey60",
        label_pal
      ),
      breaks = label_genes,       # legend shows ONLY labeled genes
      name   = "Labeled genes"
    ) +
    scale_size_identity() +
    scale_alpha_identity() +
    scale_shape_discrete(name = "Cell type") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title  = element_text(size = 20),
      axis.text   = element_text(size = 18)
    ) +
    labs(
      x = "Mean Expression (log CPM)",
      y = "Scaled Co-expression Rank"
    )
  return(p)
}


# Zinc-finger proteins
examples1 <- c(
  "ZNF384", "ZNF557", "ZNF514", "ZNF282", "ZNF333",
  "ZNF778", "ZNF383", "ZNF561", "ZNF143", "ZNF317",
  "ZNF805", "ZNF623", "ZNF583", "ZSCAN21", "ZSCAN12",
  "ZBTB17", "ZBTB26", "ZNF140", "ZNF445", "ZNF250") 

# Phosphorylation-regulated TFs
examples2 <- c(
  "ATF1", "SP2", "E2F6", 
  "CREB1", 
  #"CREM", "SP1", "STAT5B", "ATF6", "STAT1", "STAT2", "ATF2", 
  #"MEF2D" , "E2F4", "FOXP4", "ELK1", 
  "ATF6B",
   "ELK4",
   "SRF")

# Epigenetic regulators
examples3 <- c(
  "SETDB2", "MBD2", "PIAS2",
  "PIAS4", "ZBTB4", "PIAS3", "MBD1", "SETDB1")

# Hormone/signal-responsive TFs
examples4 = 
  c(
    "NR1H2","NR1H3","NR2C1","NR2C2","ARNT",
    "NRL","TP53","NR2C2",
    "GABPB1","CARF"  # less clear
  )

set.seed(123)  # for reproducibility
example_lists <- list(
  zinc_finger   = examples1,
  phosphorylation = examples2,
  epigenetic    = examples3,
  hormone       = examples4
)

plots <- lapply(names(example_lists), function(cat) {
  genes <- example_lists[[cat]]
  selected_genes <- sample(genes, min(11, length(genes)))
  
  message("Plotting category: ", cat, " | selected: ", paste(selected_genes, collapse = ", "))
  
  plot_tf_rank_vs_expr(
    tf_rpr = tf_rpr,
    tf_xpr = tf_xpr,
    label_genes = selected_genes
  ) +
    ggtitle(cat)
})

# 700 by 700 
plots[[1]]
plots[[2]]
plots[[3]]
plots[[4]]


### xpr vs rpr correlation ==============================================================
tf_xpr_avg <- tf_xpr[, .(
  mean_expr_log_cpm = mean(avg_expr_log_cpm, na.rm = TRUE)
), by = .(geneA, cell_type)] %>% rename(CellType = cell_type)

# merge
tf_merged <- tf_long %>% #tf_long is from above 
  left_join(tf_xpr_avg, by = c("geneA", "CellType"))

cor_global <- cor(tf_merged$Rank, tf_merged$mean_expr_log_cpm, method = "spearman") #spearman is better for ranks 
cor_by_cell <- tf_merged %>%
  group_by(CellType) %>%
  summarise(
    spearman_corr = cor(Rank, mean_expr_log_cpm, method = "spearman")
  )

### visualization 
plot_tf_correlation <- function(df, rank_col = "Rank", expr_col = "mean_expr_log_cpm", 
                                group_col = NULL, plot_type = c("bar", "scatter", "scatter_facet")) {
  
  plot_type <- match.arg(plot_type)
  
  rank_sym <- rlang::sym(rank_col)
  expr_sym <- rlang::sym(expr_col)
  
  if (plot_type == "bar") {
    if (is.null(group_col)) stop("group_col must be specified for bar plot")
    group_sym <- rlang::sym(group_col)
    
    cor_summary <- df %>%
      group_by(!!group_sym) %>%
      summarise(spearman_corr = cor(!!rank_sym, !!expr_sym, method = "spearman")) %>%
      ungroup() %>%
      mutate(!!group_sym := forcats::fct_reorder(!!group_sym, spearman_corr))
    
    p <- ggplot(cor_summary, aes(x = !!group_sym, y = spearman_corr)) +
      geom_col(fill = "grey30") +
      coord_flip() +
      labs(
        x = group_col,
        y = "Spearman Correlation",
        title = paste("Correlation between", expr_col, "and", rank_col)
      ) +
      my_custom_theme()
    
    return(list(summary = cor_summary, plot = p))
    
  } else if (plot_type == "scatter") {
    cor_val <- cor(df[[rank_col]], df[[expr_col]], method = "spearman")
    
    p <- ggplot(df, aes(x = !!expr_sym, y = !!rank_sym)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", color = "red") +
      labs(
        x = expr_col,
        y = rank_col,
        title = paste("Spearman correlation:", round(cor_val, 2))
      ) +
      my_custom_theme()
    
    return(list(global_correlation = cor_val, plot = p))
    
  } else if (plot_type == "scatter_facet") {
    if (is.null(group_col)) stop("group_col must be specified for faceted scatter")
    group_sym <- rlang::sym(group_col)
    
    p <- ggplot(df, aes(x = !!expr_sym, y = !!rank_sym)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", color = "red") +
      facet_wrap(vars(!!group_sym)) +
      labs(
        x = expr_col,
        y = rank_col,
        title = paste("TF Expression vs Co-expression Rank per", group_col)
      ) +
      my_custom_theme()
    
    return(list(plot = p))
  }
}

my_custom_theme <- function() {
  theme_minimal(base_size = 14) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20)
    )
}
plot_tf_correlation(tf_merged, rank_col = "Rank", expr_col = "mean_expr_log_cpm", group_col = "CellType", plot_type = "scatter_facet") #900 by 600 
plot_tf_correlation(tf_merged, rank_col = "Rank", expr_col = "mean_expr_log_cpm", group_col = NULL, plot_type = "scatter") #500 by 600 
plot_tf_correlation(tf_merged, rank_col = "Rank", expr_col = "mean_expr_log_cpm", group_col = "CellType", plot_type = "bar") #500 by 400







### variance  plots ==============================================================

#### 0. preparing the data -------------------------------------------------------
# 1) Mean expression per geneA per cell type across studies
expr_celltype_mean <- tf_xpr %>%
  group_by(geneA, cell_type) %>%
  summarise(mean_expr_ct = mean(avg_expr_log_cpm, na.rm = TRUE), .groups = "drop")

expr_stats <- expr_celltype_mean %>%
  group_by(geneA) %>%
  summarise(
    mean_expression = mean(mean_expr_ct, na.rm = TRUE),
    var_expression  = var(mean_expr_ct,  na.rm = TRUE),
    .groups = "drop"
  )

rpr_stats <- tf_rpr %>%
  group_by(geneA) %>%
  summarise(
    mean_rpr = mean(mean_ranked_rpr_scaled, na.rm = TRUE),
    var_rpr  = var(mean_ranked_rpr_scaled,  na.rm = TRUE),
    max_rpr = max(mean_ranked_rpr_scaled,  na.rm = TRUE), 
    .groups = "drop"
  )

tf_stats <- expr_stats %>% full_join(rpr_stats, by = "geneA")

TF_DBD = read.csv("/home/nelazzabi/rewiring/data/TFs/families/TF_DBD.csv")
TF_DBD = TF_DBD[, c("geneA", "DBD", "Binding.mode")]

tf_hierarchy = read.csv("/home/nelazzabi/rewiring/data/TFs/families/tf_hierarchy.csv") %>% filter(TF_isoform == "TF") %>% dplyr::select(-c("TF_isoform", "tf_name"))
tf_hierarchy$geneA = tf_hierarchy$TF_name

tf_stats <- tf_stats %>%
  left_join(TF_DBD, by = "geneA") %>%
  mutate(
    DBD = ifelse(is.na(DBD), "other", DBD),
    Binding.mode = ifelse(is.na(Binding.mode), "other", Binding.mode)
  )

tf_stats <- tf_stats %>%
  left_join(tf_hierarchy, by = "geneA") %>%
  mutate(
    superclass = ifelse(is.na(superclass), "other", superclass),
    class = ifelse(is.na(class), "other", class), 
    family = ifelse(is.na(family), "other", family),
    subfamily = ifelse(is.na(subfamily), "other", subfamily),
  )


#### 1. visualization variance -------------------------------------------------------
library(ggplot2)
library(scales)

# Discretize mean_rpr for point shape
tf_stats <- tf_stats %>%
  mutate(
    rpr_group = ntile(mean_rpr, 3),           # 3 categories: 1,2,3
    rpr_group = factor(rpr_group,
                       labels = c("low RPR", "medium RPR", "high RPR"))
  )

ggplot(tf_stats, aes(x = var_expression,
                     y = var_rpr,
                     color = mean_expression,
                     shape = rpr_group)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_c(option = "plasma", name = "Mean expression") +
  scale_shape_manual(values = c(16, 17, 15), name = "RPR group") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Variance of expression across cell types",
    y = "Variance of coexpression reproducibility across cell types"
  )
# 750 by 700




library(ggplot2)
library(ggrepel)  # for nice non-overlapping labels

plot_tf_variance <- function(tf_stats, highlight_genes = NULL) {
  
  # Separate highlighted and non-highlighted genes
  tf_stats$highlight <- tf_stats$geneA %in% highlight_genes
  
  ggplot(tf_stats, aes(x = var_expression, y = var_rpr)) +
    # all points in black first
    geom_point(data = subset(tf_stats, !highlight), 
               color = "black", shape = 16, size = 3, alpha = 0.9) +
    # highlighted points in color on top
    geom_point(data = subset(tf_stats, highlight), 
               aes(color = geneA), shape = 16, size = 3.5) +
    # labels for highlighted genes
    geom_text_repel(data = subset(tf_stats, highlight),
                    aes(label = geneA, color = geneA),
                    size = 4, show.legend = FALSE, fontface = "bold") +
    theme_minimal(base_size = 14) +
    labs(
      x = "Variance of expression across cell types",
      y = "Variance of coexpression reproducibility across cell types"
    ) +
    scale_color_manual(values = scales::hue_pal()(length(highlight_genes)))
}

plot_tf_variance(tf_stats, highlight_genes = c(""))
plot_tf_variance(tf_stats, highlight_genes = c("OLIG1", "NR2C1", "MYSM1", "OLIG2", "ZFP36L1", "ZFP64", "FIZ1", "CREM", 
                                               "RUNX1", "MEF2C", "POU2F1", "MEF2A", "FOXP1", "HLF", 
                                               "NFIB", "NFIA", "RORA", "TRPS1", 
                                               "FOS", "TSC22D4", "RFX2", "REL", 
                                               "FOXP2", "SOX6", "ZNF385B", "ZNF536"))

# 900 by 700






#### 2. categorize TFs -------------------------------------------------------
# statistically cluster TFs into regulatory regimes
library(cluster)
# Use only variance + mean properties for clustering
mat <- tf_stats %>%
  select(mean_expression, var_expression, mean_rpr, var_rpr) %>%
  scale() %>%
  as.matrix()
# Choose clusters (3–6 are common for TF regulatory states)
set.seed(1)
k <- 5  # adjust based on silhouette score if needed
cl <- kmeans(mat, centers = k, nstart = 50)
tf_stats$cluster <- factor(cl$cluster, labels = paste0("Cluster_", seq_len(k)))




# Ensure cluster is a factor -----------------------------
tf_stats$cluster <- factor(tf_stats$cluster)

# 1) Density plot — mean_expression
p1 <- ggplot(tf_stats, aes(x = mean_expression, fill = cluster, color = cluster)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  labs(title = "Distribution of mean expression by cluster",
       x = "Mean expression across cell types", y = "Density") +
  theme_minimal(base_size = 14)

# 2) Density plot — var_expression
p2 <- ggplot(tf_stats, aes(x = var_expression, fill = cluster, color = cluster)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  labs(title = "Distribution of variance in expression by cluster",
       x = "Variance of expression across cell types", y = "Density") +
  theme_minimal(base_size = 14)

# 3) Density plot — mean_rpr
p3 <- ggplot(tf_stats, aes(x = mean_rpr, fill = cluster, color = cluster)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  labs(title = "Distribution of mean reproducibility by cluster",
       x = "Mean co-expression reproducibility", y = "Density") +
  theme_minimal(base_size = 14)

# 4) Density plot — var_rpr
p4 <- ggplot(tf_stats, aes(x = var_rpr, fill = cluster, color = cluster)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  labs(title = "Distribution of variance in reproducibility by cluster",
       x = "Variance of co-expression reproducibility", y = "Density") +
  theme_minimal(base_size = 14)

library(patchwork)
combined_density <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "Distributions of TF features by regime")
# 1400 by 800 


library(uwot)
# Create matrix of standardized features -------------------
mat <- tf_stats %>%
  select(mean_expression, var_expression, mean_rpr, var_rpr) %>%
  scale() %>%
  as.matrix()
set.seed(123)
umap_out <- umap(mat, n_neighbors = 30, min_dist = 0.25, metric = "euclidean", init = "spectral")

# Add UMAP coordinates to df
tf_stats$UMAP1 <- umap_out[, 1]
tf_stats$UMAP2 <- umap_out[, 2]

p_umap <- ggplot(tf_stats, aes(x = UMAP1, y = UMAP2,
                               color = cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP of TF regulatory states in 4-D feature space",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Regulatory regime"
  )
p_umap
p_umap2 <- ggplot(tf_stats,
                  aes(x = UMAP1, y = UMAP2,
                      color = category,
                      size = mean_expression,
                      #shape = cluster
                      )) +
  geom_point(alpha = 0.85) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP of TF regulatory space with functional annotations",
    x = "UMAP 1", y = "UMAP 2",
    color = "TF category",
    size = "Mean expression",
    shape = "Regulatory regime"
  )
p_umap3 <- ggplot(tf_stats,
                  aes(x = UMAP1, y = UMAP2,
                      color = category,
                      #size = mean_expression,
                      shape = cluster
                  )) +
  geom_point(alpha = 0.85) +
  theme_minimal(base_size = 14) +
  labs(
    title = "UMAP of TF regulatory space with functional annotations",
    x = "UMAP 1", y = "UMAP 2",
    color = "TF category",
    size = "Mean expression",
    shape = "Regulatory regime"
  )
(p_umap | p_umap2) 
# 1400 by 800 





#### 3. Annotate clusters with TF biological categories -------------------------------------------------------
tf_stats <- tf_stats %>%
  mutate(category = case_when(
    geneA %in% Zinc_finger_proteins      ~ "Zinc-finger",
    geneA %in% Epigenetic_regulators ~ "Epigenetic regulator",
    geneA %in% Phosphorylation_regulated_TFs  ~ "Phosphorylation-regulated TF",
    geneA %in% Hormone_signal_responsive_TFs   ~ "Hormone/signal-responsive TF",
    TRUE ~ "Other / Unclassified"
  ))

ggplot(tf_stats, aes(x = var_expression,
                     y = var_rpr,
                     color = category,
                     shape = cluster)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Variance of expression",
    y = "Variance of reproducibility",
    color = "TF category",
    shape = "Regulatory regime",
    title = "Regulatory regimes of TFs across brain cell types"
  )


#### 4. umpas -------------------------------------------------------

library(ggplot2)
library(patchwork)
library(dplyr)
library(purrr)

# mapping of variable → legend title
color_vars <- list(
  mean_expression = "Mean expression",
  var_expression  = "Expression variance",
  mean_rpr        = "Mean coexpression rpr ranking",
  var_rpr         = "coexpression rpr ranking variance", 
  max_rpr         = "max coexpression rpr ranking", 
  DBD             = "TF DBD family", 
  Binding.mode    = "TF binding mode",
  superclass      = "TF superclass", 
  class           = "TF class", 
  family          = "TF family", 
  subfamily       = "TF subfamily"
  
)

make_plot <- function(var) {
  
  p <- ggplot(tf_stats, aes(x = UMAP1, y = UMAP2, color = .data[[var]])) +
    geom_point(alpha = 0.85) +
    theme_minimal(base_size = 14) +
    labs(
      title = color_vars[[var]],
      x = "UMAP 1", y = "UMAP 2",
      color = color_vars[[var]]
    )
  
  # choose color scale based on variable type
  if (is.numeric(tf_stats[[var]])) {
    p <- p + scale_color_viridis_c(option = "inferno")
  } else {
    p <- p + scale_color_viridis_d(option = "inferno")
  }
  
  return(p)
}


plots <- map(names(color_vars), make_plot)

# arrange (p1 | p2) / (p3 | p4)
(plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
plots[[5]]
(plots[[6]] | plots[[7]])
(plots[[8]] | plots[[9]])
(plots[[10]] | plots[[11]])
# 1400 by 800 



make_plot <- function(var, highlight = NULL, highlight_color = "red", grey_color = "grey80") {
  
  df <- tf_stats
  
  # If highlight is specified and the variable is categorical
  if (!is.null(highlight) && !is.numeric(df[[var]])) {
    df <- df %>%
      mutate(plot_color = ifelse(.data[[var]] %in% highlight, .data[[var]], "Other"))
    color_scale <- scale_color_manual(values = c(setNames(rep(highlight_color, length(highlight)), highlight), "Other" = grey_color))
    color_mapping <- "plot_color"
  } else {
    df <- df %>% mutate(plot_color = .data[[var]])
    # Automatic scale based on type
    color_scale <- if (is.numeric(df[[var]])) scale_color_viridis_c(option = "inferno") else scale_color_viridis_d(option = "inferno")
    color_mapping <- "plot_color"
  }
  
  ggplot(df, aes(x = UMAP1, y = UMAP2, color = .data[[color_mapping]])) +
    geom_point(alpha = 0.85) +
    theme_minimal(base_size = 14) +
    labs(
      title = color_vars[[var]],
      x = "UMAP 1", y = "UMAP 2",
      color = color_vars[[var]]
    ) +
    color_scale
}

make_plot("superclass", highlight = c("Zinc-coordinating DNA-binding domains"))
make_plot("class", highlight = c("Nuclear receptors with C4 zinc fingers"))

make_plot("DBD", highlight = c("C2H2 ZF"))
make_plot("Binding.mode", highlight = c("Monomer or homomultimer"))
make_plot("Binding.mode", highlight = c("Not a DNA binding protein"))












#### this is a draft for visualizing binding evidence vs [co-expression] 
#### for a given TF-gene pair relative to all other genes for that TF 



lit_targets = readRDS("/home/nelazzabi/rewiring/data/modalities/lit_curated_targets.rds")


filt_mic = mic %>% filter(geneA == "RUNX1", geneB %in% (lit_targets %>% filter(TF_Symbol == "RUNX1"))$Target_Symbol )
plot(filt_mic$avg_rank_norm_corr, filt_mic$binding_score)


mic = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Mic_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
ast = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Ast_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
opc = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Opc_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
oli = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Oli_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
exc = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Exc_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
inh = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Inh_tf_allTrg_avgRnkCoexpr_binding_lit.rds")


library(ggplot2)
library(ggrepel)
library(dplyr)

plot_geneA_binding_multi <- function(df_list, gene_pair, df_names = NULL) {
  geneA_query <- gene_pair[1]
  geneB_highlight <- gene_pair[2]
  
  if (is.null(df_names)) {
    df_names <- names(df_list)
  }
  
  combined_df <- bind_rows(
    lapply(seq_along(df_list), function(i) {
      df <- df_list[[i]]
      df$CellType <- df_names[i]
      df
    })
  )
  
  df_sub <- combined_df %>%
    filter(geneA == geneA_query) %>%
    mutate(highlight = ifelse(geneB == geneB_highlight, "highlight", "other"))
  
  df_other <- df_sub %>% filter(highlight == "other")
  df_highlight <- df_sub %>% filter(highlight == "highlight")
  
  ggplot() +
    geom_point(data = df_other, aes(x = avg_rank_norm_corr, y = binding_score),
               color = "black", alpha = 0.4, size = 1.5) +
    geom_point(data = df_highlight, aes(x = avg_rank_norm_corr, y = binding_score),
               color = "firebrick", size = 2.5) +
    geom_text_repel(data = df_highlight,
                    aes(x = avg_rank_norm_corr, y = binding_score, label = geneB),
                    color = "firebrick", fontface = "bold", size = 3.5, max.overlaps = Inf) +
    facet_wrap(~CellType, scales = "fixed") +
    labs(
      title = paste("Binding vs Coexpression for", geneA_query),
      x = "Average Normalized Rank Correlation",
      y = "Binding Score"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 13, face = "bold"),
      axis.text.y = element_text(size = 13, face = "bold")
    )
}



plot_geneA_binding_multi(
  df_list = list(inh = inh, exc = exc, ast = ast, mic = mic, opc = opc, oli = oli),
  gene_pair = c("NFIB", "HDAC4")
)

plot_geneA_binding_multi(
  df_list = list(inh = inh, exc = exc, ast = ast, mic = mic, opc = opc, oli = oli),
  gene_pair = c("MEF2C", "ANKRD44")
)

plot_geneA_binding_multi(
  df_list = list(inh = inh, exc = exc, ast = ast, mic = mic, opc = opc, oli = oli),
  gene_pair = c("JUN", "ENC1")
)























##### AUROC and AUPRC =============================================================================================================
#install.packages("precrec")

library(precrec)

# Convert labels to binary
summary_df <- exc %>% filter(geneA == "SP1") %>% 
  mutate(label = ifelse(TF_target_predicted == "positive", 1, 0))

# Check label distribution (must have both 0 and 1)
table(summary_df$label)

# Create a 'mmdata' object with scores and labels
mm_obj <- mmdata(scores = summary_df$avg_rank_norm_corr,
                 labels = summary_df$label,
                 modnames = "Aggregate Coexpression")

# Compute AUROC and AUPRC
eval_mod <- evalmod(mm_obj)

# Extract AUROC and AUPRC values
perf_summary <- as.data.frame(auc(eval_mod))
print(perf_summary)


# For coexpression
coexpr_scores <- evalmod(scores = summary_df$avg_rank_norm_corr,
                         labels = summary_df$label)

# For binding (remove NA values)
bind_df <- summary_df %>% filter(!is.na(binding_score))
bind_scores <- evalmod(scores = bind_df$binding_score,
                       labels = bind_df$label)


# Plot ROC and PR curves for coexpression
autoplot(coexpr_scores) + 
  ggtitle("Coexpression (ROC + PR)")

# Plot ROC and PR curves for binding
autoplot(bind_scores) + 
  ggtitle("Binding (ROC + PR)")











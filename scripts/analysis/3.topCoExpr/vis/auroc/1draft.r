

library(dplyr)
library(precrec)
library(purrr)
library(ggplot2)
library(patchwork)
library(ggplot2)

run_eval <- function(df, tf) {
  df_sub <- df %>%
    filter(geneA == tf) %>%
    mutate(label = ifelse(TF_target_predicted == "positive", 1, 0))
  
  # must contain both positive + negative targets
  if (length(unique(df_sub$label)) < 2) return(NULL)
  
  # evaluate ROC + PR for coexpression
  coexpr_eval <- evalmod(
    scores = df_sub$avg_rank_norm_corr,
    labels = df_sub$label
  )
  
  # evaluate ROC + PR for binding
  bind_df <- df_sub %>% filter(!is.na(binding_score))
  bind_eval <- evalmod(
    scores = bind_df$binding_score,
    labels = bind_df$label
  )
  
  list(coexpr = coexpr_eval, binding = bind_eval)
}

### apply function 
tf <- "SP1"   # <-- choose any TF
results <- map(df_list, ~ run_eval(.x, tf))
names(results) <- names(df_list)

# extract performance metrics 
performance <- imap_dfr(results, ~ {
  if (is.null(.x)) return(NULL)
  auc_df <- as.data.frame(auc(.x$coexpr))
  auc_df$CellType <- .y
  auc_df
})
performance


### plot 
plots <- imap(results, ~ {
  if (is.null(.x)) return(NULL)
  autoplot(.x$coexpr) + ggtitle(paste("Coexpression —", .y))
})
wrap_plots(plots)


plots_bind <- imap(results, ~ {
  if (is.null(.x)) return(NULL)
  autoplot(.x$binding) + ggtitle(paste("Binding —", .y))
})
wrap_plots(plots, plots_bind, ncol = 2)















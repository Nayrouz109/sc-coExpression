


library(precrec)
library(purrr)
library(dplyr)

tf <- "SP1"
tf <- "NFIB"

scores_list <- list()
labels_list <- list()
model_names <- c()


### coexpression 
for (cell in names(df_list)) {
  df <- df_list[[cell]] %>% 
    filter(geneA == tf) %>% 
    mutate(label = ifelse(TF_target_predicted == "positive", 1, 0))
  
  # skip datasets lacking both 0 and 1 labels
  if (length(unique(df$label)) < 2) next
  
  scores_list[[cell]] <- df$avg_rank_norm_corr
  labels_list[[cell]] <- df$label
  model_names <- c(model_names, cell)
}

mm_combined <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = model_names
)

eval_combined <- evalmod(mm_combined)
autoplot(eval_combined) + ggtitle(paste("Coexpression ROC + PR —", tf))


### binding 
scores_list <- list()
labels_list <- list()
model_names <- c()

for (cell in names(df_list)) {
  df <- df_list[[cell]] %>% 
    filter(geneA == tf, !is.na(binding_score)) %>% 
    mutate(label = ifelse(TF_target_predicted == "positive", 1, 0))
  
  if (length(unique(df$label)) < 2) next
  
  scores_list[[cell]] <- df$binding_score
  labels_list[[cell]] <- df$label
  model_names <- c(model_names, cell)
}

mm_combined_bind <- mmdata(
  scores = scores_list,
  labels = labels_list,
  modnames = model_names
)

eval_combined_bind <- evalmod(mm_combined_bind)
autoplot(eval_combined_bind) + ggtitle(paste("Binding ROC + PR —", tf))




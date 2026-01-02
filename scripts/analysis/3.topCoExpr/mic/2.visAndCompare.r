
###  pick TFs  ===========================================================================================
dir <- "/home/nelazzabi/rewiring/data/coExpr/AD-Ctl"
files <- list.files(dir, pattern = ".*tf_reprod_tbl.rds", full.names = TRUE)

data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))
data_list[["Mic_tf_reprod_tbl_rnk"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/rnkAgg", "Mic_tf_reprod_tbl.rds")) 

data_list[["Mic_tf_reprod_tbl_sfl"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/shuffled", "Mic_tf_reprod_tbl.rds")) 
data_list[["Mic_tf_reprod_tbl_sfl_tr"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/shuffled/tr", "Mic_tf_reprod_tbl.rds")) 
data_list[["Mic_tf_reprod_tbl_dsp"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/downsampling/AD-Ctl-smpls", "Mic_tf_reprod_tbl.rds"))
data_list[["Mic_tf_reprod_tbl_dsc"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/downsampling/AD-Ctl-smpls-clls", "Mic_tf_reprod_tbl.rds"))
data_list[["Mic_tf_reprod_tbl_wm"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/wm", "Mic_tf_reprod_tbl.rds"))
data_list[["Mic_tf_reprod_tbl_dsp_ws"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/downsampling/AD-Ctl-smpls-ws", "Mic_tf_reprod_tbl.rds"))
data_list[["Mic_tf_reprod_tbl_dsc_ws"]] = readRDS(file.path("/home/nelazzabi/rewiring/data/coExpr/downsampling/AD-Ctl-smpls-clls-ws", "Mic_tf_reprod_tbl.rds"))

tf_reprod_avg = readRDS("/home/nelazzabi/rewiring/data/coExpr/downsampling/AD-Ctl-dsp-n/Mic_tf_reprod_tbl.rds")

tf_reprod_combined <- bind_rows(tf_reprod_avg, .id = "iterations")
data_list[["Mic_tf_reprod_tbl_avg"]] <- tf_reprod_combined %>%
  group_by(geneA) %>%
  summarise(
    avg_meanTop200_AD = mean(meanTop200_AD, na.rm = TRUE),
    avg_meanTop200Ctl = mean(meanTop200Ctl, na.rm = TRUE),
    mic_TF = any(mic_TF),  # If geneA is classified as mic_TF in any dataset, keep TRUE
    brain_TF = any(brain_TF)  # If geneA is classified as brain_TF in any dataset, keep TRUE
  ) %>%
  ungroup()

data_list$Mic_tf_reprod_tbl_avg$meanTop200_AD = data_list$Mic_tf_reprod_tbl_avg$avg_meanTop200_AD
data_list$Mic_tf_reprod_tbl_avg$meanTop200Ctl = data_list$Mic_tf_reprod_tbl_avg$avg_meanTop200Ctl 

### visualize TF reproducibaility ***************************************************
int_genes = (data_list$Mic_tf_reprod_tbl %>% mutate(delta = abs(meanTop200Ctl - meanTop200_AD)) %>% filter(delta >= 10))$geneA
plot_tf_reprod(data_list$Ast_tf_reprod_tbl, int_genes)
plot_tf_reprod(data_list$Exc_tf_reprod_tbl, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_rnk, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_avg, int_genes)

plot_tf_reprod(data_list$Mic_tf_reprod_tbl_sfl, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_sfl_tr, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_dsp, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_dsc, int_genes)

plot_tf_reprod(data_list$Mic_tf_reprod_tbl_wm, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_dsp_ws, int_genes)
plot_tf_reprod(data_list$Mic_tf_reprod_tbl_dsc_ws, int_genes)


plot_tf_reprod(tf_reprod, int_genes)

### visualize top reproducible TF ***************************************************
plot_tf_summary0(data_list$Mic_tf_reprod_tbl, "mic_TF", "brain_TF") 
plot_tf_summary0(data_list$Mic_tf_reprod_tbl_dsp, "mic_TF", "brain_TF") 
plot_tf_summary0(data_list$Mic_tf_reprod_tbl_dsc, "mic_TF", "brain_TF") 
plot_tf_summary0(data_list$Mic_tf_reprod_tbl_rnk, "mic_TF", "brain_TF")
plot_tf_summary(data_list$Mic_tf_reprod_tbl, "brain_TF") 
plot_tf_summary(data_list$Ast_tf_reprod_tbl, "brain_TF")
plot_tf_summary(data_list$Exc_tf_reprod_tbl, "brain_TF")

### focus on specific TFs ***************************************************
int_genes_top = c("RUNX2", "NFATC2", "MAF", "RUNX1", "MEF2A", "MEF2C", "E2F3", "SMAD3", "E2F3", "YBX1", "FOXN3")
int_genes_top = c("MEF2C" , "RUNX1", "HIF1A" )
plot_tf_reprod2(data_list$Mic_tf_reprod_tbl, int_genes_top, int_genes_top)
plot_tf_reprod2(data_list$Exc_tf_reprod_tbl, int_genes_top, int_genes_top)
plot_tf_reprod2(data_list$Ast_tf_reprod_tbl, int_genes_top, int_genes_top)
plot_tf_reprod2(data_list$Mic_tf_reprod_tbl_sfl, int_genes_top, int_genes_top)
plot_tf_reprod2(data_list$Mic_tf_reprod_tbl_sfl_tr, int_genes_top, int_genes_top)






###  loso  ===========================================================================================

loso = readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_tf_reprod_tbl_loso.rds")

ggplot(loso, aes(x = meanTop200_AD_LOSO, y = meanTop200_AD)) +
  geom_point(alpha = 0.4, size = 1.5) +  
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +
  facet_wrap(~ Study_Removed) +
  labs(x = "Mean Top 200 AD (LOSO)", y = "Mean Top 200 AD") +
  theme_minimal()

ggplot(loso, aes(x = meanTop200_CTL_LOSO, y = meanTop200_CTL)) +
  geom_point(alpha = 0.4, size = 3) +  
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +
  facet_wrap(~ Study_Removed) +
  labs(x = "meanTop200_CTL_LOSO", y = "meanTop200_CTL") +
  theme_minimal()

ggplot(loso, aes(x = meanTop200_AD_LOSO, y = meanTop200_CTL_LOSO)) +
  geom_point(alpha = 0.4, size = 2) +  
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "solid", linewidth = 1) +
  facet_wrap(~ Study_Removed) +
  labs(x = "Mean Top 200 AD (LOSO)", y = "Mean Top 200 CTL (LOSO)") +
  theme_minimal()


test$delta_AD <- test$meanTop200_AD - test$meanTop200_AD_LOSO
test$delta_CTL <- test$meanTop200_CTL - test$meanTop200_CTL_LOSO


ggplot(test, aes(x = delta_CTL, y = delta_AD, color = Study_Removed)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
  labs(x = "Delta CTL (meanTop200_CTL - meanTop200_CTL_LOSO)", 
       y = "Delta AD (meanTop200_AD - meanTop200_AD_LOSO)", 
       title = "Comparison of Delta AD vs Delta CTL per Study") +
  theme_minimal() +
  facet_wrap(~Study_Removed)
theme(legend.position = "bottom")



test$delta <- abs (test$meanTop200_AD_LOSO - test$meanTop200_AD)
ggplot(test, aes(x = geneA, y = delta, group = Study_Removed, color = Study_Removed)) +
  geom_point(size = 3) +
  geom_line() +
  labs(x = "Gene", y = "Delta (meanTop200_AD_LOSO - meanTop200_AD)", title = "Delta per Gene across Studies") +
  theme_minimal() +
  theme(legend.position = "bottom")





###  pick TFs - targets  ===========================================================================================

selected_tfs = (data_list$Mic_tf_reprod_tbl %>% mutate(delta = abs(meanTop200Ctl - meanTop200_AD)) %>% filter(delta >= 20))$geneA
selected_tfs <- c("RUNX1", "NFATC2", "MAF", "MEF2A", "MEF2C", "E2F3")  

mic_tf_list <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_tf_top200_cor_intg_list.rds")
merged_all_mic <- do.call(rbind, lapply(mic_tf_list[selected_tfs], \(df) df[setdiff(names(df), "SCENIC_TF2G_rho")]))


merged_all_mic_cor <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_tf_allTrg_list.rds")

TF_target_list <- merged_all_mic %>% filter(delta_abs >= 5) %>%
  group_by(geneA) %>%
  summarize(geneB_list = list(geneB), .groups = "drop") %>%
  deframe()
extract_TF_targets <- function(df, TF_target_list) {
  df %>%
    filter(geneA %in% names(TF_target_list)) %>%
    filter(mapply(`%in%`, geneB, TF_target_list[geneA])) 
}

all_TF_targets_df <- purrr::map_dfr(merged_all_mic_cor, extract_TF_targets, TF_target_list = TF_target_list, .id = "dataset") %>% 
  mutate(study = gsub("(_[^_]+)$", "", dataset))

t_test_results <- all_TF_targets_df %>%
  group_by(geneA, geneB) %>%
  filter(n_distinct(condition) > 1) %>%  # Ensure both conditions exist
  summarize(
    p_value = t.test(rank_norm_corr ~ condition)$p.value,
    mean_AD = mean(rank_norm_corr[condition == "AD"], na.rm = TRUE),
    mean_CTL = mean(rank_norm_corr[condition == "control"], na.rm = TRUE),
    diff_mean = mean_AD - mean_CTL,
    n_AD = sum(condition == "AD"),
    n_CTL = sum(condition == "control"),
    .groups = "drop"
  ) %>%
  arrange(p_value)


merged_all_ast_cor <- readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Ast_tf_allTrg_list.rds")
all_TF_targets_df_ast <- purrr::map_dfr(merged_all_ast_cor, extract_TF_targets, TF_target_list = TF_target_list, .id = "dataset") %>% 
  mutate(study = gsub("(_[^_]+)$", "", dataset))

plot_rank_norm_corr(all_TF_targets_df %>% filter(geneA == "MAF"), "NINJ1")
plot_rank_norm_corr(all_TF_targets_df_ast %>% filter(geneA == "MAF"), "NINJ1")


plot_rank_norm_corr(all_TF_targets_df %>% filter(geneA == "NFATC2"), "FBXO11")
plot_rank_norm_corr(all_TF_targets_df_ast %>% filter(geneA == "NFATC2"), "FBXO11")


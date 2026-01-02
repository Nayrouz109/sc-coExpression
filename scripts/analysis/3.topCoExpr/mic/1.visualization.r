


##### read in data ========================================================================================================================
##### =====================================================================================================================================


#### 1. TFs 
TFs = readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS") 
TFs <- lapply(TFs, function(df) {df %>% mutate(geneB = Symbol)})
lit_TFs = fread("/home/nelazzabi/rewiring/data/TFs/Curated_targets_all_Sept2023.tsv")
#alex_mic_allRank  = readRDS("/home/nelazzabi/rewiring/data/TFs/alex_mic_allRank.rds")
#alex_mic_fishZ    = readRDS("/home/nelazzabi/rewiring/data/TFs/alex_mic_fishZ.rds")


eRegulon = fread("/home/nelazzabi/rewiring/data/TFs/eRegulon_metadata_filtered.csv")
eRegulonV1 = eRegulon %>% filter(TF %in% gene_list$transcription_factors) %>% 
  dplyr::select(c("TF", "Gene")) %>%  unique() 

freq = as.data.frame(table(eRegulonV1$Gene))
length(unique((freq %>% filter(Freq > 1))$Var1))



microglia_TFs_chatGPT <- c("SPI1", "IRF8", "SALL1", "MEF2C", "MafB")
microglia_TFs_paper <- colnames(test)[-1]
mic_TF = unique(c(microglia_TFs_chatGPT, microglia_TFs_paper))




#### 2. AD-genes  
mic_neuroinflammation = fread("/home/nelazzabi/rewiring/data/AD-genes/AD_genes_mic_neuroinflammation.csv")
SEA_AD_micGenes = fread("/home/nelazzabi/rewiring/data/AD-genes/SEAAD_mic_genes.csv") 



#### 2. my data 
my_data_TF_top200 = readRDS("/cosmos/data/project-data/NW-rewiring/coExpr/disease/TF-targets/tf_top200_list.rds")
my_data_TF_allTrg = readRDS("/cosmos/data/project-data/NW-rewiring/coExpr/disease/TF-targets/tf_allTrg_list.rds") 




##### creat freq table ========================================================================================================================
##### =====================================================================================================================================
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/freqTbls.r")
source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/processTFsTargets.R")


freqTable = summarize_gene_pairs(my_data_TF_top200)
tf_summary_df <- summarize_tf_targets(freqTable)
tf_reprod <- calculate_tf_target_reproducibility(my_data_TF_top200)


for (i in 1:length(names(my_data_TF_allTrg))) {
  my_data_TF_allTrg[[i]]$condition = ifelse(grepl("_AD", names(my_data_TF_allTrg)[i]), "AD" , "control" ) 
}
combined_tf_top200 <- bind_rows(my_data_TF_top200, .id = "study")
combined_tf_allTrg <- bind_rows(my_data_TF_allTrg, .id = "study")
rm(my_data_TF_allTrg)


tf_targets_alex_fishZ  = rank_gene_partners (readRDS("/home/amorin/scratch/aggregate_microglia/Cormats/Hg_pcor/aggregate_cormat_FZ_hg.RDS")[["Agg_mat"]], "RUNX1")
tf_targets_alex_allRank = rank_gene_partners (readRDS("/home/amorin/scratch/aggregate_microglia/Cormats/Hg_pcor/aggregate_cormat_allrank_filter_lenient_hg.RDS")[["Agg_mat"]], "RUNX1")




##### visualization ========================================================================================================================
##### =====================================================================================================================================

# 1. bar plots freq distribution 
ggplot(freqTable %>% filter(freq_CTL > 0), aes(x = freq_CTL)) +
  geom_bar(fill = "blue", alpha = 0.8) +
  labs(title = "Distribution of TF-top 200 co-expr Pairs Across Studies",
       x = "Number of Studies",
       y = "Number of TF-target Pairs") +
  scale_x_continuous(breaks = 1:14) +  # Ensure x-axis has integer labels from 1 to 14
  theme_classic()


# 2.1 dot plots - TF reproducability 
ggplot(tf_summary_df, aes(x = proportion_reproducible_CTL, y = reorder(geneA, proportion_reproducible_CTL))) +
  geom_point(color = "black", size = 2) +  # Dot plot with black points
  labs(
    x = "Proportion of geneB targets reproducible across studies",
    y = "TF (geneA)"
  ) +
  theme_classic() + 
  theme(
    axis.text.y = element_blank())
#So, if proportion_reproducible = 0.6, then 60% of the targets for that TF are observed 
# in at least two studies, making them reproducible 
# 20% of MEF2A’s geneB targets appear in 6 or more studies. 

# 2.2 dot plots - TF unique targets 
ggplot(tf_summary_df, aes(x = reorder(geneA, total_targets_CTL), y = total_targets_CTL)) +
  geom_point(color = "black", size = 2) +  # Dot plot
  labs(
    x = "TF (geneA)",
    y = "Frequency (Number of geneB targets)") + 
  theme_classic() + 
  theme(axis.text.x = element_blank())



# 2.3  dot plots TF reproducability (fancy version)
library(ggrepel)  
highlight_tfs <- c("FOXN3", "MEF2A", "RUNX1", "FOXP2", "JAZF1")

ggplot(tf_summary_df, aes(x = proportion_reproducible_CTL, y = reorder(geneA, proportion_reproducible_CTL))) +
  geom_point(aes(color = geneA %in% highlight_tfs), size = 2) +  # Color specific TFs
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +  # Red for selected TFs, black for others
  geom_text_repel(
    data = tf_summary_df %>% filter(geneA %in% highlight_tfs),  # Only label the highlighted TFs
    aes(label = geneA),
    color = "red",
    size = 4,  # Adjust label size
    nudge_x = 0.02,  # Small horizontal shift to avoid overlap
    direction = "both"
  ) +
  labs(
    x = "Proportion of geneB targets reproducible across studies",
    y = "TF (geneA)"
  ) +
  theme_classic()


# 3. correlation distributions CTL vs AD
ggplot(RUNX1_cor, aes(x = rank_norm_corr, color = condition)) +
  geom_density(alpha = 0.5) + 
  scale_color_manual(values = c("black", "firebrick")) + 
  facet_wrap(~study, scales = "free") +  # Separate plots per study
  theme_minimal() +
  labs(title = "Distribution of Correlation Values per Study (Matching TF-Target Pairs)",
       x = "Rank Normalized Correlation",
       y = "Density",
       fill = "Condition") 


# 4.1 scatter plot nairuz vs Alex's - FZ 
ggplot(RUNX1_cor %>% filter(condition == "control"), aes(x = rank_norm_corr, y = coExpr_rank_FZ)) +
  geom_point(alpha = 0.1, color = "blue", size = 1.5) +
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "Comparison of my data vs. alex's (FZ) Rank Correlations - controls only",
    x = "my data",
    y = "alex's data"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))

# 4.2 scatter plot nairuz vs Alex's - allRanks 
#AD
#control
ggplot(RUNX1_cor %>% filter(condition == "control"), aes(x = rank_norm_corr, y = coExpr_rank_allRank)) +
  geom_point(alpha = 0.1, color = "blue", size = 1.5) +
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "Comparison of my data vs. alex's (allRank) Rank Correlations - controls only",
    x = "my data",
    y = "alex's data"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))


ggplot(RUNX1_cor, aes(x = rank_norm_corr, y = coExpr_rank_FZ, color = condition)) +
  geom_point(alpha = 0.5, size = 1.5) +  # Adjust transparency and size for better visibility
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "Comparison of my data vs. Alex's (allRank) Rank Correlations",
    x = "My data",
    y = "Alex's data",
    color = "Condition"  # Legend label
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))







##### biological analysis - logic ========================================================================================================================
##### ============================================================================================================================================
# 1) intersect these with Alex’s microglial coexpression for these TFs and 
# 2) intersect with the literature-curated targets for these TFs. 
# 3) whether they are not DE in AD.



# 1 is helpful as confirmation (at least, for the ones that have lower coexp in AD than controls — 
#                              or the ones that have lower coexp in controls should not show up in Alex’s data as well?)
# 2 is helpful to make it biologically relevant
# 3. to evaluate if DE’s a likely confound.

#The most interesting would check all three boxes obviously. 
# I notice that NFATC2 has RUNX1 as a partner on your list, and vice versa. That is, they are both TRs and they are showing differential coexpression with each other pretty strongly, 
# stronger in AD. Something like that might be worth looking at a little more closely, especially as both are mentioned in the paper.



##### biological analysis - RUNX1 \ NFATC2 \ MAF ========================================================================================================================
##### ============================================================================================================================================

NFATC2_merged_all = integrate_gene_ranking(freqTable, 
                                           tf_targets_alex_fishZ, 
                                           tf_targets_alex_allRank, 
                                           TFs, 
                                           "NFATC2", 
                                           lit_TFs, 
                                           NULL, 
                                           mic_neuroinflammation, 
                                           SEA_AD_micGenes, 
                                           eRegulon) 

MAF_merged_all = integrate_gene_ranking(freqTable, 
                                           tf_targets_alex_fishZ, 
                                           tf_targets_alex_allRank, 
                                           TFs, 
                                           "MAF", 
                                           lit_TFs, 
                                           NULL, 
                                           mic_neuroinflammation, 
                                           SEA_AD_micGenes, 
                                           eRegulon) 
MAF_rnk_integrated_cor_all <- integrate_freq_ranking(gene_name = "MAF", combined_tf_allTrg = combined_tf_allTrg, RUNX1_merged_all = MAF_merged_all)






RUNX1_merged_all = integrate_gene_ranking(freqTable, 
                                          tf_targets_alex_fishZ, 
                                          tf_targets_alex_allRank, 
                                          TFs, 
                                          "RUNX1", 
                                          lit_TFs, 
                                          NULL, 
                                          mic_neuroinflammation, 
                                          SEA_AD_micGenes, 
                                          eRegulon) 

RUNX2_merged_all = integrate_gene_ranking(freqTable, 
                                          tf_targets_alex_fishZ, 
                                          tf_targets_alex_allRank, 
                                          TFs, 
                                          "RUNX2", 
                                          lit_TFs, 
                                          NULL, 
                                          mic_neuroinflammation, 
                                          SEA_AD_micGenes, 
                                          eRegulon) 


E2F3_merged_all = integrate_gene_ranking(freqTable, 
                                          tf_targets_alex_fishZ, 
                                          tf_targets_alex_allRank, 
                                          TFs, 
                                          "E2F3", 
                                          lit_TFs, 
                                          NULL, 
                                          mic_neuroinflammation, 
                                          SEA_AD_micGenes, 
                                          eRegulon) 
E2F3_rnk_integrated_cor_all <- integrate_freq_ranking(gene_name = "E2F3", combined_tf_allTrg = combined_tf_allTrg, RUNX1_merged_all = E2F3_merged_all)


MEF2C_merged_all = integrate_gene_ranking(freqTable, 
                                          tf_targets_alex_fishZ, 
                                          tf_targets_alex_allRank, 
                                          TFs, 
                                          "MEF2C", 
                                          lit_TFs, 
                                          NULL, 
                                          mic_neuroinflammation, 
                                          SEA_AD_micGenes, 
                                          eRegulon) 








RUNX1_rnk_integrated_cor <- filter_and_integrate_ranking(
  delta_threshold = 5,  # Set to NULL if no filtering is needed
  gene_name = "RUNX1",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = RUNX1_merged_all
)

NFATC2_rnk_integrated_cor <- filter_and_integrate_ranking(
  delta_threshold = 5,  
  gene_name = "NFATC2",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = NFATC2_merged_all
)

MAF_rnk_integrated_cor <- filter_and_integrate_ranking(
  delta_threshold = 5,  
  gene_name = "MAF",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = MAF_merged_all
)

ggplot(MAF_rnk_integrated_cor, aes(x = log10(rank_norm_corr + 1), y = log10(coExpr_rank_allRank +1), color = condition)) +
  geom_point(alpha = 0.5, size = 1.5) +  # Adjust transparency and size for better visibility
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "MAF comparison of my data vs. Alex's (allRank) Rank Correlations",
    x = "My data",
    y = "Alex's data",
    color = "Condition"  # Legend label
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))





RUNX1_rnk_integrated_cor_top <- filter_and_integrate_ranking(
  ctl_threshold = 5,  # Set to NULL if no filtering is needed
  gene_name = "RUNX1",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = RUNX1_merged_all
)

NFATC2_rnk_integrated_cor_top <- filter_and_integrate_ranking(
  ctl_threshold = 5, 
  gene_name = "NFATC2",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = NFATC2_merged_all
)

MAF_rnk_integrated_cor_top <- filter_and_integrate_ranking(
  ctl_threshold = 5, 
  gene_name = "MAF",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = MAF_merged_all
)

ggplot(RUNX1_rnk_integrated_cor_top, aes(x = log10(rank_norm_corr + 1), y = log10(coExpr_rank_allRank +1), color = condition)) +
  geom_point(alpha = 0.5, size = 1.5) +  # Adjust transparency and size for better visibility
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "RUNX1 comparison of my data (top CTL) vs. Alex's (allRank) Rank Correlations",
    x = "My data",
    y = "Alex's data",
    color = "Condition"  # Legend label
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))




### extracting correlatons 
RUNX1_rnk_integrated_cor_all <- integrate_freq_ranking(
  gene_name = "RUNX1",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = RUNX1_merged_all
)

RUNX2_rnk_integrated_cor_all <- integrate_freq_ranking(
  gene_name = "RUNX2",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = RUNX2_merged_all
)

NFATC2_rnk_integrated_cor_all <- integrate_freq_ranking(
  gene_name = "NFATC2",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = NFATC2_merged_all
)

MEF2C_rnk_integrated_cor_all <- integrate_freq_ranking(
  gene_name = "MEF2C",
  combined_tf_allTrg = combined_tf_allTrg,
  RUNX1_merged_all = MEF2C_merged_all
)





####### interesting genes for the NIH grant 
genes <- c("FOXN3", "FCHSD2", "RUNX1", 
           "FBXO11", "MERTK")
plots <- lapply(genes, function(g) plot_rank_norm_corr(NFATC2_rnk_integrated_cor_all, g))

genes <- c("BHLHE41", "STX7")
plots <- lapply(genes, function(g) plot_rank_norm_corr(MAF_rnk_integrated_cor_all, g))

genes <- c("ITPR2", "TRPM2")
plots <- lapply(genes, function(g) plot_rank_norm_corr(E2F3_rnk_integrated_cor_all, g))

genes <- c("ELL2", "GLS", 
           "NFATC2", "SRGAP2C", "TANC1")
plots <- lapply(genes, function(g) plot_rank_norm_corr(RUNX1_rnk_integrated_cor_all, g))

wrap_plots(plots,  nrow = 1)  # Arrange in a 2x2 grid








ggplot(MAF_rnk_integrated_cor_all %>% filter(condition == "control"), aes(x = rank_norm_corr, y = freq_CTL)) +
  geom_point(alpha = 0.3, size = 2) +  
  facet_wrap(~study, scales = "free") + 
  theme_minimal() +
  labs(
    title = "MAF",
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))

RUNX1_diff <- calculate_differences(RUNX1_rnk_integrated_cor_all)
NFATC2_diff <- calculate_differences(NFATC2_rnk_integrated_cor_all)
MAF_diff <- calculate_differences(MAF_rnk_integrated_cor_all)

plot_delta_relationship(RUNX1_diff, title = "RUNX1")
plot_delta_relationship(RUNX1_diff, title = "RUNX1", highlight_geneB = "PIK3CD")
plot_delta_relationship(NFATC2_diff, title = "NFATC2")
plot_delta_relationship(MAF_diff, title = "MAF")





































##### TFs with top delta  ========================================================================================================================

top_1_percent <- data_list$Mic_tf_reprod_tbl %>%
  mutate(delta = abs(meanTop200Ctl - meanTop200_AD)) %>%
  arrange(desc(delta)) %>%
  slice_head(prop = 0.10) %>%  # Select top 1%
  select(geneA)
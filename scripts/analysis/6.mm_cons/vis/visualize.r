

### visualization of mouse cons. signal 

allRank = readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_CTL_geneA_hg_mm_sCor_summary.rds") 
FZ      = readRDS("/home/nelazzabi/rewiring/data/coExpr/AD-Ctl/Mic_CTL_mm_FZ_geneA_hg_mm_sCor_summary.rds")





ggplot(allRank, aes(x = top_200, y = s.cor)) +
  geom_point(aes(color = ifelse(is_TF, "firebrick", "black"), 
                 size = ifelse(is_TF, 3.5, 2),
                 alpha = ifelse(is_TF, 0.6, 0.2)),  # Moved inside aes()
             shape = 1, stroke = 1) +  # Shape 1 makes it an open circle
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    x = "GeneA",
    y = "Spearman Correlation",
    title = "Correlation between Nairuz and Alex Rank Normalized Correlations"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_classic() +
  ylab("Spearman's cor") + 
  xlab("Top200 in CTL") + 
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        legend.position = "none") +
  scale_color_identity() +  # Uses raw colors in aes(color = ...)
  scale_size_identity() +   # Uses raw size values
  scale_alpha_identity()    # Ensures alpha is mapped correctly


library(ggplot2)
library(ggrepel)


files <- list.files(dir, pattern = ".*rpr.*.rds", full.names = TRUE)
data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))
top_10 = (tail(data_list$Mic_t200_rpr_tf %>% arrange(meanTop200Ctl), 20))$geneA

ggplot(allRank, aes(x = top_200, y = s.cor)) +
  geom_point(aes(color = ifelse(is_TF, "firebrick", "black"), 
                 size = ifelse(is_TF, 2, 1), 
                 alpha = ifelse(is_TF, 0.6, 0.3)), 
             shape = 1, stroke = 1) +  # Open circles
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  
  # Add text labels for selected genes with lines
  geom_text_repel(data = subset(allRank, geneA %in% top_10),
                  aes(label = geneA), size = 5, box.padding = 0.5, max.overlaps = 30) +
  
  theme_minimal() +
  labs(
    x = "GeneA",
    y = "Spearman Correlation",
    title = "Correlation between Nairuz and Alex Rank Normalized Correlations"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_classic() +
  ylab("Spearman's cor") + 
  xlab("Top200 in CTL") + 
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        legend.position = "none") +
  scale_color_identity() +  # Uses raw colors in aes(color = ...)
  scale_size_identity()     # Uses raw size values




# Assign colors
allRank$color <- with(allRank, ifelse(is_TF, "firebrick",
                                      ifelse(is_ribo, "orange", "black")))

# Assign size: emphasize TF or ribo
allRank$size <- ifelse(allRank$is_TF | allRank$is_ribo, 2, 1)

# Assign alpha: emphasize TF or ribo
allRank$alpha <- ifelse(allRank$is_TF | allRank$is_ribo, 0.6, 0.3)


ggplot(allRank, aes(x = top_200, y = s.cor)) +
  geom_point(aes(color = color, 
                 size = size, 
                 alpha = ifelse(is_TF | is_ribo, 0.8, 0.1)
                 ), 
             shape = 1, stroke = 1) +  # Open circles
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  
  geom_text_repel(data = subset(allRank, geneA %in% top_10),
                  aes(label = geneA), size = 5, box.padding = 0.5, max.overlaps = 30) +
  
  theme_minimal() +
  labs(
    x = "GeneA",
    y = "Spearman Correlation",
    title = "Correlation between Nairuz and Alex Rank Normalized Correlations"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_classic() +
  ylab("Spearman's cor") + 
  xlab("Top200 in CTL") + 
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20),
        legend.position = "none") +
  scale_color_identity() +
  scale_size_identity() +
  scale_alpha_identity()

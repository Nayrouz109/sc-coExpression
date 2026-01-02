

#### disease ========================================================
all_dge = readRDS("/cosmos/data/project-data/NW-rewiring/dge/AD2/all_res.rds")

microglia_dge_df <- bind_rows(
  lapply(names(all_dge), function(study_name) {
    df <- all_dge[[study_name]]$Mic
    df <- df %>% mutate(Study = study_name)  # Add study name as a column
    return(df)
  })
)

#### cell-types ========================================================
all_dge = readRDS("/cosmos/data/project-data/NW-rewiring/data/3.prcsd-slcGenes/PB/cellType-integrated/DEA/all_res.rds")$Mic_Exc

for (i in 1:length(all_dge)) {
  print(dim(all_dge[[i]]))
  
}

microglia_dge_df <- bind_rows(
  lapply(names(all_dge), function(study_name) {
    df <- all_dge[[study_name]]
    df <- df %>% mutate(Study = study_name)  # Add study name as a column
    return(df)
  })
)


ggplot(microglia_dge_df, aes(x = P.Value)) +
  #geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_histogram(binwidth = 0.03, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~ Study, scales = "free_y") +  # Facet by Study with independent y-axes
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +  # Force x-axis from 0 to 1
  theme_classic() +
  labs(title = "P-value Distribution Across Studies",
       x = "P-value",
       y = "Frequency") +
  theme(strip.text = element_text(size = 10, face = "bold"))



threshold_line <- -log10(0.05)
ggplot(microglia_dge_df, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.2, color = "black") +  
  geom_hline(yintercept = threshold_line, color = "firebrick", size = 0.5) +  # Threshold line
  facet_wrap(~ Study, scales = "free") +  # Facet by Study
  theme_classic() +
  labs(title = "Faceted Volcano Plot of Microglia DGE",
       x = "logFC",
       y = "-log10(Adjusted P-value)") +
  theme(strip.text = element_text(size = 10, face = "bold"))




### this is an addition to TF_vis.r script 
### this is to help prioritize TFs 


library(dplyr)
library(purrr)
library(stringr)

# DATA --------------------------------------------------------------------------------------------------------------------------------------------------
dir <- "/home/nelazzabi/rewiring/data/coExpr/0all-human/1.1TFtop200-coExpr-perCell"
files <- list.files(dir, pattern = ".*rpr.*.rds", full.names = TRUE)

data_list <- setNames(lapply(files, readRDS), tools::file_path_sans_ext(basename(files)))
data_list <- data_list[grep("rpr_tf", names(data_list))]

tf_long <- imap_dfr(data_list, ~
                      mutate(.x,
                             CellType = str_remove(.y, "_t200_rpr_tf"))
) %>% select(-c("meanTop200_AD"))

# Weighted reproducibility table per TF × CellType -----------------------------------------------------------------------------------------------------
tf_weighted <- tf_long %>%
  group_by(CellType) %>%
  mutate(
    Rank = rank(-meanTop200Ctl, ties.method = "min"),
    #WeightedScore = max(Rank) - Rank + 1,   # converts rank to a weight (higher rank = higher weight)
    #NormScore = scales::rescale(meanTop200Ctl, to = c(0, 1))  # optional normalization
  ) %>%
  ungroup()

### create threshold table 
tf_thresholds <- tf_weighted %>%
  mutate(
    top30  = Rank <= 30,
    top50  = Rank <= 50,
    top100 = Rank <= 100,
    top150 = Rank <= 150,
    top200 = Rank <= 200,
    WeightedScore = (top30 * 5) + (top50 * 4) + (top100 * 3) + (top150 * 2) + (top200 * 1)
  ) %>%
  select(geneA, CellType, starts_with("top"), WeightedScore) %>% 
  mutate(across(starts_with("top"), as.integer))

# normalize weighted score per cell type 
tf_summary <- tf_thresholds %>%
  group_by(CellType) %>%
  mutate(
    NormScore = WeightedScore / max(WeightedScore, na.rm = TRUE)
  ) %>%
  ungroup()


# global summary 
global_summary <- tf_summary %>%
  mutate(in_top50 = top50 == 1) %>%
  group_by(geneA) %>%
  summarise(
    `#CellTypes_in_top50` = sum(in_top50, na.rm = TRUE),
    MeanWeightedScore = mean(WeightedScore, na.rm = TRUE),
    cell_type = paste(unique(CellType[in_top50]), collapse = ", ")
  ) %>%
  mutate(
    Category = case_when(
      `#CellTypes_in_top50` >= 5 ~ "Global core",
      #`#CellTypes_in_top50` == 1 ~ paste0(cell_type, "-specific"),
      `#CellTypes_in_top50` == 1 ~ "Cell-type specific",
      TRUE ~ "Shared regulators"
    )
  ) %>%
  arrange(desc(MeanWeightedScore)) %>%
  filter(`#CellTypes_in_top50` > 0)













### visualization ------------------------------------------------------------------------------------------------------------------------------

# 1 350 by 2000 
### take the union of the top50 for all 6 cell types 
### change the plots such that y - axis is TF and x-axis is the cell type and arrange plots by cell type 
library(dplyr)
library(ggplot2)
library(viridis)

tf_union_scores <- tf_summary %>%
  filter(top200 > 0) %>%
  group_by(CellType) %>%
  arrange(desc(NormScore)) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  distinct(geneA) %>%
  inner_join(tf_summary, by = "geneA") %>%
  mutate(geneA = fct_reorder(geneA, NormScore, mean, .desc = TRUE))

ggplot(tf_union_scores, aes(x = CellType, y = geneA, fill = NormScore)) +
  geom_tile() +
  scale_fill_viridis_c(option = "cividis") +
  labs(
    x = "Cell Type",
    y = "Transcription Factor",
    fill = "Norm. Reproducibility",
    title = "Reproducibility of Top TFs Across Brain Cell Types"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))








# 2 600 by 2000
library(forcats)
ggplot(global_summary, aes(x = fct_rev(reorder(geneA, MeanWeightedScore)),
                           y = MeanWeightedScore, fill = Category)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Global core" = "#FDE725",
                               "Shared regulators" = "#A79D74",
                               "Cell-type specific" = "#081F4A")) +
  labs(x = "Transcription Factor", y = "Mean Weighted Score",
       title = "Classification of TFs by reproducibility patterns") +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 15))


# 400 by 600 
ggplot(global_summary, aes(x = `#CellTypes_in_top50`, fill = Category)) +
  geom_histogram(binwidth = 1, color = "white") +
  scale_fill_manual(values = c(
    "Global core" = "#FDE725",
    "Shared regulators" = "#A79D74",
    "Cell-type specific" = "#081F4A"
  )) +
  labs(
    x = "Number of cell types (TF in top50)",
    y = "Count of TFs"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 15),
    legend.position = "none"   # removes the legend
  )






# 3 500 by 600 
### vis MeanWeightedScore per TF category 
library(ggplot2)

ggplot(global_summary, aes(x = MeanWeightedScore, fill = Category)) +
  geom_density(alpha = 0.8, color = NA) +
  scale_fill_manual(values = c(
    "Global core" = "#FDE725",
    "Shared regulators" = "#A79D74",
    "Cell-type specific" = "#081F4A"
  )) +
  labs(
    #title = "Distribution of TF Mean Weighted Reproducibility Scores",
    x = "Mean Weighted Score",
    y = "Density",
    fill = "TF Category"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_text(size = 15),
    legend.position = "top"
  )


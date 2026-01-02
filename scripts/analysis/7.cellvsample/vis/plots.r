
library(dplyr)
library(ggplot2)
library(patchwork)

####  functions ================================================================================================
####  functions ================================================================================================
####  functions ================================================================================================
plot_reproducibility <- function(df, reprod_label) {
  # 1. filter the dataframe for the given reproducibility label
  df_filtered <- df %>%
    filter(reprod == reprod_label)
  
  # 2. compute mean odds ratios
  line_df <- df_filtered %>%
    group_by(cell_type, reference_dataset) %>%
    summarize(mean_or = mean(odds_ratio, na.rm = TRUE), .groups = "drop")
  
  # 3. create the plot
  plot <- ggplot(df_filtered, aes(x = reference_dataset, y = odds_ratio, color = cell_type
                                  )) +
    # per–cell type points (hollow)
    geom_point(shape = 21, fill = NA, size = 2) +
    
    # mean line across datasets (only this connects)
    geom_line(data = line_df,
              aes(x = reference_dataset, y = mean_or, group = cell_type),
              linewidth = 1) +
    
    # big filled mean points on top
    geom_point(data = line_df,
               aes(x = reference_dataset, y = mean_or, color = cell_type
                   ),
               size = 4) +
    
    theme_minimal(base_size = 12) +
    labs(
      x = "Dataset",
      y = "OR",
      title = paste("Reproducibility Across Datasets:", reprod_label)
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_y_continuous(trans = "log2") #+
    #scale_y_continuous(trans = "log2", limits = c(1, 4050))
  
  return(plot)
}
plot_reproducibility <- function(df, reprod_label) {
  df_filtered <- df %>%
    filter(reprod == reprod_label)
  
  line_df <- df_filtered %>%
    group_by(cell_type, reference_dataset) %>%
    summarize(mean_or = mean(odds_ratio, na.rm = TRUE), .groups = "drop")
  
  plot <- ggplot(df_filtered, aes(x = reference_dataset, y = odds_ratio)) +
    # per–dataset dots, filled by dataset
    geom_point(aes(fill = dataset), shape = 21, color = "black", size = 3) +
    
    # mean line across datasets (colored by cell type)
    geom_line(data = line_df,
              aes(x = reference_dataset, y = mean_or, group = cell_type, color = cell_type),
              linewidth = 1) +
    
    # big filled mean points on top (colored by cell type)
    geom_point(data = line_df,
               aes(x = reference_dataset, y = mean_or, color = cell_type),
               size = 4) +
    
    theme_minimal(base_size = 12) +
    labs(
      x = "Dataset",
      y = "OR",
      title = paste("Reproducibility Across Datasets:", reprod_label)
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 12)
    ) +
    scale_y_continuous(trans = "log2")
  
  return(plot)
}
plot_comparison <- function(df, comparison) {
  stopifnot(length(comparison) == 2)
  
  comp1 <- comparison[1]
  comp2 <- comparison[2]
  
  # compute mean OR per group
  df_means <- df %>%
    group_by(cell_type, reference_dataset, reprod) %>%
    summarize(mean_or = mean(odds_ratio, na.rm = TRUE), .groups = "drop")
  
  # split into the two groups
  df1 <- df_means %>%
    filter(reprod == comp1) %>%
    rename(mean_or_1 = mean_or)
  
  df2 <- df_means %>%
    filter(reprod == comp2) %>%
    rename(mean_or_2 = mean_or)
  
  # join on reference_dataset + cell_type
  df_compare <- inner_join(
    df1,
    df2,
    by = c("reference_dataset", "cell_type")
  )
  
  # build labels dynamically
  xlab_txt <- paste0(comp1, " co-expression reproducibility\n(Mean odds ratio - log2 scale)")
  ylab_txt <- paste0(comp2, " co-expression reproducibility\n(Mean odds ratio - log2 scale)")
  title_txt <- paste0("Reproducibility of ", comp2, " Co-expression\nvs. ", comp1, " Co-expression")
  
  # plot
  plot <- ggplot(df_compare, aes(x = mean_or_1, y = mean_or_2)) +
    geom_point(aes(color = reference_dataset, shape = cell_type),
               show.legend = TRUE, size = 6) +
    geom_abline(slope = 1, linetype = "dashed", color = "gray40") +
    scale_x_continuous(trans = "log2", limits = c(1, 600)) +
    scale_y_continuous(trans = "log2", limits = c(1, 600)) +
    labs(
      x = xlab_txt,
      y = ylab_txt,
      title = title_txt
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 20),
      strip.text = element_text(size = 14),
      plot.title = element_text(size = 23),
      axis.title = element_text(size = 18)
    )
  
  return(plot)
}

##### DATA ================================================================================================
xCell      = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/0all/rpr/fisher/xCellxSubj_allCell_types_repro_fisher_CTL.csv") %>% mutate(reprod = "xCell")
xCellxSubj = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher/xCellxSubj_allCell_types_repro_fisher_CTL.csv") %>% mutate(reprod = "xCellxSubj")
xSubj      = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xSubject/rpr/fisher/xCellxSubj_allCell_types_repro_fisher_CTL.csv") %>% mutate(reprod = "xSubj")
xSubjxType = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType/rpr/fisher/xCellxSubj_allCell_types_repro_fisher_CTL.csv")  %>% mutate(reprod = "xSubjectxType")


xCell_all      = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/0all/rpr/fisher/xCellxSubj_allCell_types_repro_fisherN_ALL.csv") %>% mutate(reprod = "xCell")
xCellxSubj_all = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher/xCellxSubj_allCell_types_repro_fisherN_ALL.csv") %>% mutate(reprod = "xCellxSubj")
xSubj_all      = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xSubject/rpr/fisher/xCellxSubj_allCell_types_repro_fisherN_ALL.csv") %>% mutate(reprod = "xSubj")
#xSubjxType_all = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType/rpr/fisher/xCellxSubj_allCell_types_repro_fisher_CTL.csv")  %>% mutate(reprod = "xSubjectxType")



bxu_all <- read_csv('/space/scratch/bxu_project/pavlidislab_project/pan/final/2-reprod/odd_ratios.csv')


xSubjxTypeN = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xSubjectxType/rpr/fisher/xCellxSubj_allCell_types_repro_fisherN_CTL.csv")  %>% mutate(reprod = "xSubjectxType")  
xCellxSubjects1 = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher/xCellxSubj_all_cell_types_repro_fisher.csv")
xCellxSubjects2 = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher/xCellxSubj_all_cell_types_repro_fisher_CTL.csv")
xCell_nairuz = read.csv("/cosmos/data/project-data/NW-rewiring/coExpr/0all/rpr/fisher/xCellxSubj_all_cell_types_repro_fisher_CTL.csv")

combined_df <- df_all %>%
  select(reference_dataset, dataset, cell_type, odds_ratio, pvalue, reprod) %>%
  bind_rows(xCellxSubjects)



##### ==================================================================================================================================================================
##### brianns ======================================================================================================================================
##### ==================================================================================================================================================================

bxu_cell <- plot_reproducibility(bxu_all, "xcell")
bxu_subj <- plot_reproducibility(bxu_all, "xsubj")
y_range <- range(bxu_all$odds_ratio, na.rm = TRUE) # Find global y-range across both plots
bxu_cell <- bxu_cell + coord_cartesian(ylim = y_range) # Apply same limits
bxu_subj <- bxu_subj + coord_cartesian(ylim = y_range) # Apply same limits
bxu_cell + bxu_subj + plot_layout(ncol = 2, guides = "collect")


##### ==================================================================================================================================================================
##### nairuz's work [23 datasets] ======================================================================================================================================
##### ==================================================================================================================================================================

p1 <- plot_reproducibility(xCell,       "xCell")
p2 <- plot_reproducibility(xCellxSubj,  "xCellxSubj")
p3 <- plot_reproducibility(xSubj,       "xSubj")
p4 <- plot_reproducibility(xSubjxType,  "xSubjectxType")

y_range_log2 <- range((c(
    xCell$odds_ratio,
    xCellxSubj$odds_ratio,
    xSubj$odds_ratio,
    xSubjxType$odds_ratio
  )),
  na.rm = TRUE
)

p1 <- p1 + coord_cartesian(ylim = y_range_log2)
p2 <- p2 + coord_cartesian(ylim = y_range_log2)
p3 <- p3 + coord_cartesian(ylim = y_range_log2)
p4 <- p4 + coord_cartesian(ylim = y_range_log2)

(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

p1 + p2 + plot_layout(ncol = 2, guides = "collect")
p3 + p4 + plot_layout(ncol = 2, guides = "collect")



##### ==================================================================================================================================================================
##### nairuz's work [18 datasets] ======================================================================================================================================
##### ==================================================================================================================================================================
datasets_main <- c(
  "ling-2024", "ROSMAP-Mathys-2024", "MSBB-2024", "RADC-2024", "HBCC-2024", "SZBDMulti-2024",
  "MSSM-2024", "Rexach-2024", "ROSMAP-erosion-2023", "Kamath-NPH-2023",
  "Hoffman-2023", "ROSMAP-Mathys-2019", "ROSMAP-Mathys-2023", "ROSMAP-TREM2-2020",
  "ROSMAP-DeJager-2024", "ROSMAP-DeJager-2023", "SEA-AD-MTG-2024", "SEA-AD-DLPFC-2024"
)

datasets_other <- c(
  "Lau-2020", "Lim-2022", "Nagy-2020", "Pineda-2021", "Velmeshev-2019"
)

p11 <- plot_reproducibility( xCell %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main), "xCell")
p22 <- plot_reproducibility( xCellxSubj %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main), "xCellxSubj") 
p33 <- plot_reproducibility( xSubj %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main), "xSubj")
p44 <- plot_reproducibility( xSubjxType %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main), "xSubjectxType")

y_range_log2 <- range((c(
    (xCell       %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main))$odds_ratio,
    (xCellxSubj  %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main))$odds_ratio,
    (xSubj       %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main))$odds_ratio,
    (xSubjxType  %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main))$odds_ratio
  )),
  na.rm = TRUE
)

p1 <- p1 + coord_cartesian(ylim = y_range_log2)
p2 <- p2 + coord_cartesian(ylim = y_range_log2)
p3 <- p3 + coord_cartesian(ylim = y_range_log2)
p4 <- p4 + coord_cartesian(ylim = y_range_log2)

p1 + p2 + plot_layout(ncol = 2, guides = "collect")
p3 + p4 + plot_layout(ncol = 2, guides = "collect")




##### ==================================================================================================================================================================
##### global y-range ======================================================================================================================================
##### ==================================================================================================================================================================

datasets_list <- list(xCell, xCellxSubj, xSubj, xSubjxType, bxu_all)
y_range_global <- range(unlist(lapply(datasets_list, function(df) {(df$odds_ratio)})), na.rm = TRUE)

bxu_cell <- bxu_cell + coord_cartesian(ylim = y_range_global) # Apply same limits
bxu_subj <- bxu_subj + coord_cartesian(ylim = y_range_global) # Apply same limits
p1 <- p1 + coord_cartesian(ylim = y_range_global)
p2 <- p2 + coord_cartesian(ylim = y_range_global)
p3 <- p3 + coord_cartesian(ylim = y_range_global)
p4 <- p4 + coord_cartesian(ylim = y_range_global)

p11 <- p11 + coord_cartesian(ylim = y_range_global)
p22 <- p22 + coord_cartesian(ylim = y_range_global)
p33 <- p33 + coord_cartesian(ylim = y_range_global)
p44 <- p44 + coord_cartesian(ylim = y_range_global)

bxu_cell + bxu_subj + plot_layout(ncol = 2, guides = "collect")
p1 + p2 + plot_layout(ncol = 2, guides = "collect")
p3 + p4 + plot_layout(ncol = 2, guides = "collect")
p11 + p22 + plot_layout(ncol = 2, guides = "collect")
p33 + p44 + plot_layout(ncol = 2, guides = "collect")



plot_reproducibility(xCellxSubjects1, "xcellxSubj")
plot_reproducibility(xCellxSubjects2, "xcellxSubj")
plot_reproducibility(combined_df, "xcellxSubj")
plot_reproducibility(combined_df, "xcell")
plot_reproducibility(combined_df, "xsubj")
plot_reproducibility(xCell_nairuz, "xcellxSubj")



########################################## XCELL vs XSUBJ REPROD #############################################################################

combined_df <- bind_rows(xCell, xCellxSubj, xSubj, xSubjxType)
glimpse(combined_df)

plot_comparison(combined_df %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main) , c("xCell", "xCellxSubj"))
plot_comparison(combined_df %>% filter(reference_dataset %in% datasets_main & dataset %in% datasets_main) , c("xCell", "xSubj"))
plot_comparison(combined_df, c("xCell", "xCellxSubj"))
plot_comparison(combined_df, c("xCell", "xSubj"))

plot_comparison(bxu_all, c("xcell", "xsubj"))

plot_comparison(combined_df, c("xcell", "xsubj"))
plot_comparison(combined_df, c("xcell", "xcellxSubj"))
plot_comparison(combined_df, c("xsubj", "xcellxSubj"))














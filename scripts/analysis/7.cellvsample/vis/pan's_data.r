


##### ==================================================================================================================================================================
##### pan's data [CTL vs ALL] ======================================================================================================================================
##### ==================================================================================================================================================================
datasets_other <- c(
  "Lau-2020", "Lim-2022", "Nagy-2020", "Pineda-2021", "Velmeshev-2019", "ROSMAP-Mathys-2019"
)

p1all <- plot_reproducibility( xCell_all , "xCell")
p2all <- plot_reproducibility( xCellxSubj_all , "xCellxSubj") 
p3all <- plot_reproducibility( xSubj_all , "xSubj")
p1ctl <- plot_reproducibility( xCell %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other), "xCell")
p2ctl <- plot_reproducibility( xCellxSubj %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other), "xCellxSubj")
p3ctl <- plot_reproducibility( xSubj %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other), "xSubj")

y_range_pans <- range((c(
  (xCell       %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other))$odds_ratio,
  (xCellxSubj  %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other))$odds_ratio,
  (xSubj       %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other))$odds_ratio,
  (xCell_all)$odds_ratio,
  (xCellxSubj_all)$odds_ratio, 
  (xSubj_all)$odds_ratio
)),
na.rm = TRUE
)

p1all <- p1all + coord_cartesian(ylim = y_range_pans)
p2all <- p2all + coord_cartesian(ylim = y_range_pans)
p3all <- p3all + coord_cartesian(ylim = y_range_pans)
p1ctl <- p1ctl + coord_cartesian(ylim = y_range_pans)
p2ctl <- p2ctl + coord_cartesian(ylim = y_range_pans)
p3ctl <- p3ctl + coord_cartesian(ylim = y_range_pans)

p1all + p2all + plot_layout(ncol = 2, guides = "collect")
p1ctl + p2ctl + plot_layout(ncol = 2, guides = "collect")

p1all + p1ctl + plot_layout(ncol = 2, guides = "collect")
p2all + p2ctl + plot_layout(ncol = 2, guides = "collect")
p3all + p3ctl + plot_layout(ncol = 2, guides = "collect")

combined_pans_all <- bind_rows(xCell_all, 
                           xCellxSubj_all, 
                           xSubj_all) 
combined_pans_ctl <- bind_rows(xCell %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other), 
                               xCellxSubj %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other), 
                               xSubj %>% filter(reference_dataset %in% datasets_other & dataset %in% datasets_other)) 
                            
plot_comparison(combined_pans_all, c("xCell", "xCellxSubj"))
plot_comparison(combined_pans_all, c("xCell", "xSubj"))
plot_comparison(combined_pans_ctl, c("xCell", "xCellxSubj"))
plot_comparison(combined_pans_ctl, c("xCell", "xSubj"))




mm_data = fread("/cosmos/data/project-data/NW-rewiring/data/mouse/0.original/mm_study_summary.csv")


# Create the lollipop plot
ggplot(mm_data, aes(x = Num_Cells / 1000, y = reorder(Study, Num_Cells))) +
  geom_segment(aes(xend = 0, yend = Study), color = "darkgrey") +
  geom_point(color = "orange", size = 2) +
  labs(x = "Number of Cells (in thousands)", y = "Study", title = "Lollipop Plot of Cell Counts per Study") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 5))


library(ggplot2)

ggplot(mm_data, aes(x = Num_Cells, y = reorder(Study, Num_Cells))) +
  geom_segment(aes(xend = 1, yend = Study), color = "gray") +
  geom_point(color = "blue", size = 3) +
  scale_x_log10() +  # Log-scale for better visualization
  labs(x = "Number of Cells (log10)", y = "Study", title = "Lollipop Plot of Cell Counts per Study") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 5))

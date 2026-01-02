
### function ===================================================================
combine_summary <- function(dir_path) {
  # List all files matching the pattern
  files <- list.files(path = dir_path, pattern = "*_corSummary.csv", full.names = TRUE)
  
  # Read and process each file
  summary_list <- lapply(files, function(file) {
    filename <- basename(file)
    
    # Split the filename by underscore
    parts <- strsplit(filename, "_")[[1]]
    study_name <- parts[2]
    condition <- parts[length(parts) - 1]
    cell_type <- parts[1]
    
    # Read the summary file
    dt <- fread(file)
    
    # Add new columns
    dt[, study_name := study_name]
    dt[, condition := condition]
    dt[, cell_type := cell_type]
    
    # Reorder columns
    dt <- dt[, .(study_name, condition, total_genes, total_single_cells, total_patients, cell_type)]
    return(dt)
  })
  
  # Combine all summaries
  combined_summary <- rbindlist(summary_list, fill = TRUE)
  return(combined_summary)
}





### implementation ===================================================================

dsp_n_dir = "/cosmos/data/project-data/NW-rewiring/coExpr/downsample/AD-Ctl-dsp-n/cor-summary"
dir_path_dsp <- "/cosmos/data/project-data/NW-rewiring/coExpr/downsample/AD-Ctl-smpls/cor-summary"
dir_path_dsc <- "/cosmos/data/project-data/NW-rewiring/coExpr/downsample/AD-Ctl-smpls-clls/cor-summary"
dir_path_sfl_tr <- "/cosmos/data/project-data/NW-rewiring/coExpr/shuffled/tr/cor-summary"

dir_path <- "/cosmos/data/project-data/NW-rewiring/coExpr/disease/cor-summary"

combined_summary_mic = combine_summary(dir_path) %>% filter(cell_type == "Mic")
combined_summary_ast = combine_summary(dir_path) %>% filter(cell_type == "Ast")
combined_summary_exc = combine_summary(dir_path) %>% filter(cell_type == "Exc")

combined_summary_mic_dsp_n = combine_summary(dsp_n_dir) #%>% mutate(cell_type = paste0(cell_type, "_dsp"))
combined_summary_mic_dsp = combine_summary(dir_path_dsp) %>% mutate(cell_type = paste0(cell_type, "_dsp"))
combined_summary_mic_dsc = combine_summary(dir_path_dsc)
combined_summary_mic_sfl_tr = combine_summary(dir_path_sfl_tr) %>% mutate(cell_type = paste0(cell_type, "_sfl_tr"))


# Write the combined summary to a CSV file in the same directory (or change the path as needed)
output_file1 <- file.path(dir_path, "Mic_combined_summary.csv")
output_file2 <- file.path(dir_path, "Ast_combined_summary.csv")
output_file3 <- file.path(dir_path, "Exc_combined_summary.csv")
fwrite(combined_summary_mic, file = output_file1)
fwrite(combined_summary_ast, file = output_file2)
fwrite(combined_summary_exc, file = output_file3)




### visualization ===================================================================
library(dplyr)
library(ggplot2)
library(tidyr)

# Combine all datasets
combined_summary <- bind_rows(combined_summary_mic, 
                              combined_summary_ast, 
                              combined_summary_exc)

combined_summary <- bind_rows(combined_summary_mic, 
                              combined_summary_mic_dsp, 
                              combined_summary_mic_sfl_tr)

# Calculate AD-to-CTL ratio per study per cell type
summary_ratios <- combined_summary %>%
  group_by(study_name, cell_type, condition) %>%
  summarise(total_single_cells = sum(total_single_cells), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = total_single_cells) %>%
  mutate(AD_CTL_ratio = AD / CTL)

summary_ratios_patient_norm <- combined_summary %>%
  group_by(study_name, cell_type, condition) %>%
  summarise(total_single_cells = sum(total_single_cells),
            total_patients = sum(total_patients), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = c(total_single_cells, total_patients)) %>%
  mutate(avg_cells_per_patient_AD = total_single_cells_AD / total_patients_AD,
         avg_cells_per_patient_CTL = total_single_cells_CTL / total_patients_CTL,
         AD_CTL_patient_norm_ratio = avg_cells_per_patient_AD / avg_cells_per_patient_CTL)


summary_ratios <- combined_summary %>%
  group_by(study_name, cell_type, condition) %>%
  summarise(total_single_cells = sum(total_single_cells),
            total_patients = sum(total_patients), .groups = 'drop') %>%
  pivot_wider(names_from = condition, values_from = c(total_single_cells, total_patients)) %>%
  mutate(AD_CTL_patient_ratio = total_patients_AD / total_patients_CTL,
         avg_cells_per_patient_AD = total_single_cells_AD / total_patients_AD,
         avg_cells_per_patient_CTL = total_single_cells_CTL / total_patients_CTL,
         AD_CTL_patient_norm_ratio = avg_cells_per_patient_AD / avg_cells_per_patient_CTL)

### 1.  sn AD to CTL ratio 
ggplot(summary_ratios, aes(x = study_name, y = AD_CTL_ratio, group = cell_type, color = cell_type)) +
  geom_point(size = 5, shape = "diamond") +
  geom_line(linewidth = 1, 
    #linetype = "dotted", lwd = 1
    ) +
  scale_color_manual(values = c("Mic" = "blue", "Ast" = "firebrick", "Exc" = "chartreuse4")) +
  labs(title = "sn AD to CTL Ratio Trend Across Studies",
       y = "AD:CTL Ratio",
       color = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size = 10)) 

### 2.  patient AD to CTL ratio 
ggplot(summary_ratios, aes(x = study_name, y = AD_CTL_patient_ratio, group = cell_type, color = cell_type)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1, alpha = 0.8) +
  #scale_color_manual(values = c("Mic" = "blue", "Ast" = "firebrick", "Exc" = "chartreuse4")) +
  scale_color_manual(values = c("Mic" = "blue", "Mic_dsp" = "firebrick", "Mic_sfl_tr" = "chartreuse4")) +
  labs(title = "AD to CTL Patient Ratio Across Studies",
       y = "AD:CTL Patient Ratio",
       color = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size = 10)) 

### 3.  sn AD to CTL ratio - patient normalized 
ggplot(summary_ratios, aes(x = study_name, y = AD_CTL_patient_norm_ratio, group = cell_type, color = cell_type)) +
  geom_point(size = 5) +
  geom_line(linewidth = 1.5, alpha = 1) +
  #scale_color_manual(values = c("Mic" = "#000075", "Ast" = "#b0b0ff", "Exc" = "blue")) +
  scale_color_manual(values = c("Mic" = "blue", "Ast" = "#eb8100", "Exc" = "chartreuse4")) +
  #scale_color_manual(values = c("Mic" = "blue", "Mic_dsp" = "firebrick", "Mic_sfl_tr" = "chartreuse4")) +
  labs(title = "Patient-Normalized AD to CTL Single Cell Ratio Across Studies",
       x = "Study",
       y = "Patient-Normalized AD:CTL Single Cell Ratio",
       color = "Cell Type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        axis.text.y = element_text(size = 20), 
        axis.title = element_text(size = 10))






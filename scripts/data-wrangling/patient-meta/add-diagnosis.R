

library(data.table)

# Read the data
df1 <- as.data.frame(
  data.table::fread('/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData.csv')
)

# Update the 'disease' column based on 'study' and 'diagnosis'
df1$disease <- with(df1, ifelse(
  study %in% c("ling-2024", "MSSM-2024", "SZBDMulti-2024") & diagnosis == "yes", "SCZ",
  ifelse(study %in% c("ling-2024", "MSSM-2024", "SZBDMulti-2024") & diagnosis == "no", "CTL",
         ifelse(study %in% c("Kamath-NPH-2023", "ROSMAP-DeJager-2024", "ROSMAP-erosion-2023", "ROSMAP-DeJager-2023", 
                             "Rexach-2024", "Hoffman-2023", "ROSMAP-TREM2-2020", "ROSMAP-Mathys-2019", "ROSMAP-Mathys-2024", 
                             "SEA-AD-MTG-2024", "SEA-AD-DLPFC-2024", "ROSMAP-Mathys-2023") & diagnosis == "yes", "AD",
                ifelse(study %in% c("Kamath-NPH-2023", "ROSMAP-DeJager-2024", "ROSMAP-erosion-2023", "ROSMAP-DeJager-2023", 
                                    "Rexach-2024", "Hoffman-2023", "ROSMAP-TREM2-2020", "ROSMAP-Mathys-2019", "ROSMAP-Mathys-2024", 
                                    "SEA-AD-MTG-2024", "SEA-AD-DLPFC-2024", "ROSMAP-Mathys-2023") & diagnosis == "no", "CTL",
                       ifelse(study %in% c("MSBB-2024", "HBCC-2024", "RADC-2024") & disease == "control", "CTL", 
                              disease))))))

df1 <- df1 %>% select(-c(1,2)) %>% 
  mutate(disease = case_when(
    study == "SZBDMulti-2024" & grepl("BD", patientID) ~ "BD",
    study == "SZBDMulti-2024" & grepl("CON", patientID) ~ "CTL",
    study == "SZBDMulti-2024" ~ "SCZ",
    TRUE ~ disease  # Retain existing disease value for other studies
  ))

disease_counts <- df1 %>%
  distinct(patientID, study, disease, .keep_all = TRUE) %>%  # Keep unique patientID per study and disease
  group_by(study, disease) %>%
  summarise(count = n(), .groups = 'drop')


# Save the updated dataframe
output_path <- '/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv'
write.csv(df1, output_path, row.names = FALSE)





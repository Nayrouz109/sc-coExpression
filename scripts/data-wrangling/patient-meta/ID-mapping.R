

#### updating the psych patient IDs using the mapping file provided by authors ===========================
RADC-2024
HBCC-2024
MSBB-2024
patients <- data.table::fread("/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv", header = TRUE)
psych_p <- patients %>% filter(study %in% c("RADC-2024", "HBCC-2024", "MSBB-2024"))
psych_ID_mapping <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/psychAD-2024/data/meta/NPS-AD_DonorID_mapping.csv", header = TRUE)
merged_data <- merge(psych_p, psych_ID_mapping, by.x = "patientID", by.y = "DonorID", all.x = TRUE)


psych_pF = merged_data %>% 
  mutate(patientID = individualID) %>%  # Replace patientID with individualID
  select(-individualID)  # Remove the individualID column




#### updating ROSMAP-Mathys 2023 ===========================
Mathys-2019
ROSMAP-DeJager-2024

Mathys-2023

patients_ROSMAP = patients %>% filter(study == "ROSMAP-Mathys-2023")
individual_metadata = data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-patient-metaData/MIT-multiome/MIT_ROSMAP_Multiomics_individual_metadata.csv") %>% select(c("individualID", "subject"))
ROSMAP_clinical <- read_csv("/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/meta/ROSMAP_clinical.csv") %>% select(c("projid", "individualID"))


patients_ROSMAP_mapped <- patients_ROSMAP %>%
  left_join(individual_metadata, by = c("patientID" = "subject")) %>%  # Join on subject
  left_join(ROSMAP_clinical, by = "individualID") %>%  # Join on individualID for additional data
  unique()  # Remove duplicate rows, if needed

patients_ROSMAP_mapped$patientID <-  patients_ROSMAP_mapped$projid 
patients_ROSMAP_mapped <- patients_ROSMAP_mapped %>% select(-c("individualID", "projid"))


#### updating all patient meta file ===========================
patients <- data.table::fread("/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv", header = TRUE)
patients_V <- patients %>% filter(!(study %in% c("ROSMAP-Mathys-2023", "RADC-2024", "HBCC-2024", "MSBB-2024"))) 
patients_VF <- bind_rows(patients_V, patients_ROSMAP_mapped, psych_pF)  



### sanity check ============================================
# make sure study names are consistent between file names and study names in patient metadata file 
matrices_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw"
matrices <- list.files(matrices_dir, pattern = "h5ad")
study_names <-unique(sub("^[^_]+_(.*)\\.h5ad$", "\\1", matrices))

study_names[(!(study_names %in% patients_VF$study))]

write.csv(patients_VF, "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData.csv", row.names = TRUE)


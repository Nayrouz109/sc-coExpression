




### ROSMAP - MIT =========================================================================================================================== 
meta_patient <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT/data/0-source/meta/individual_metadata_deidentified.tsv")
meta_patient_final <- meta_patient %>%  mutate(sex = ifelse (msex == 0, "female", "male"), 
                                               brainRegion = "PFC", 
                                               study = "ROSMAP-MIT", 
                                               ADdiag3types = Pathologic_diagnosis_of_AD) %>% 
                                        select(c("subject", "sex", "age_death", "Pathologic_diagnosis_of_AD", "ADdiag3types" , "brainRegion", "study")) 


colnames(meta_patient_final) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
saveRDS(meta_patient_final, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT/data/0-source/meta/ROSMAP-MIT-patient.rds" )


### ROSMAP- Micro =========================================================================================================================== 
ROSMAP_micro_patient <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-Micro/data/0-source/meta/ROSMAP-MICRO.csv")

ROSMAP_micro_patient_final <- ROSMAP_micro_patient %>%  mutate(sex = ifelse (msex == "Female", "female", "male"), 
                                                               study = "ROSMAP-microGlia",
                                                               ADdiag2types = ifelse(ADdiag2types == "nonAD", "no", "yes")) %>% 
  select(c("subject", "sex", "age_death", "ADdiag2types", "ADdiag3types" , "brainRegion", "study")) 


colnames(ROSMAP_micro_patient_final) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
saveRDS(ROSMAP_micro_patient_final, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Micro/data/0-source/meta/ROSMAP-Micro-patient.rds" )


### ROSMAP- erosion =========================================================================================================================== 

patient_erosion <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-erosion/data/0-source/meta/erosion_patient.csv")

patient_erosion_final <- patient_erosion %>%  mutate(sex = ifelse (sex == "Female", "female", "male"), 
                                                     study = "ROSMAP-erosion", 
                                                     brainRegion = "PFC", 
                                                     ADdiag2types = ifelse(AD_group == "nonAD", "no", "yes")) %>% 
  select(c("subject", "sex", "age_death", "ADdiag2types", "AD_group" , "brainRegion", "study"))

colnames(patient_erosion_final) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
saveRDS(patient_erosion_final, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-erosion/data/0-source/meta/ROSMAP-Erosion-patient.rds" )




### ROSMAP- gStability =========================================================================================================================== 

gStability_p <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-Gstability/data/0-source/meta/human_smart-seq2_fusions.maxsensitivity_prediction_numbers.meta.counts.tsv")

gStability_p <- gStability_p %>% select(c("subject", "age_death", "msex", "cogdx.ad", "braak56")) %>% unique() %>% 
  mutate(sex = ifelse (msex == 0, "female", "male"), 
         study = "ROSMAP-gStability", 
         brainRegion = "PFC", 
         ADdiag2types = ifelse(braak56 == "CTRL", "no", "yes")) %>% 
  select(c("subject", "sex", "age_death", "ADdiag2types", "cogdx.ad" , "brainRegion", "study"))
  
colnames(gStability_p) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
saveRDS(gStability_p, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Gstability/data/0-source/meta/ROSMAP-gStability-patient.rds" )


### ROSMAP- DeJager ===========================================================================================================================
de_Jager_meta <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager/data/0-source/meta/cell-annotation.full-atlas.csv", header = TRUE)
clinical <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/MIT-ROSMAP/meta/ROSMAP_clinical.csv")

de_Jager_annotations_clinical  <- left_join(de_Jager_meta[, c("individualID")], clinical)
de_Jager_annotations_clinical <- de_Jager_annotations_clinical[!(duplicated(de_Jager_annotations_clinical)) ,]

### reported controls | ADs 
# non-AD = 159 
# AD     = 265 

dim(filter(de_Jager_annotations_clinical, ceradsc %in% c("4", "3") & braaksc %in% c("0", "1", "2", "3", "4")))   

deJager_p <- de_Jager_annotations_clinical %>% 
  mutate(sex = ifelse (msex == 0, "female", "male"), 
         study = "ROSMAP-DeJager",
         brainRegion = "DLPFC", 
         ADdiag2types = ifelse(ceradsc %in% c("4", "3") & braaksc %in% c("0", "1", "2", "3", "4"), "no", "yes"), 
         ADdiag3types =  ADdiag2types ) %>% 
  select(c("individualID", "sex", "age_death", "ADdiag2types", "ADdiag3types" , "brainRegion", "study"))

colnames(deJager_p) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
saveRDS(deJager_p, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager/data/0-source/meta/ROSMAP-DeJager-patient.rds" )

### ROSMAP- Mathys1 ===========================================================================================================================
mathys2019 <- readRDS("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/Mathys-2019/data/original-big-matrix/Mathys_2019_data_ready.rds")
mathys2019p <- mathys2019$meta %>% select(c("projid", "msex", "pathology.group", "AD.class")) %>% unique() %>% 
  mutate(sex = msex, 
         study = "ROSMAP-Mathys1", 
         brainRegion = "PFC",  
         ADdiag2types = ifelse(AD.class == "control", "no", "yes"), 
         ADdiag3types = pathology.group) %>% 
  select(c("projid", "sex", "ADdiag2types", "ADdiag3types" , "brainRegion", "study"))

colnames(mathys2019p) = c("patientID", "sex", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
rownames(mathys2019p) = NULL
saveRDS(mathys2019p, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/meta/ROSMAP-Mathys1-patient.rds" )



### Kamath- NPH ===========================================================================================================================
kamath_dir = "/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/Kamath/annotations/rds-list"
kamath_data <- readRDS(file.path(kamath_dir, "kamath-meta.rds"))

NPH_meta <- bind_rows(lapply(kamath_data, function(df) { 
  dplyr::filter(df, anno_organism == "human" & ds_batch == "human_NPH")
}))

NPH_pateint <- NPH_meta %>% select(c("anno_batch", "anno_sex", "anno_age",  "anno_condition", "anno_region")) %>% unique() %>% 
  mutate(sex = ifelse(anno_sex == "M", "male", "female"), 
         study = "Kamath-NPH", 
         brainRegion = anno_region,  
         ADdiag2types = ifelse(anno_condition == "Ctrl", "no", "yes"), 
         ADdiag3types = anno_condition, 
         anno_batch = gsub("human_", "", anno_batch)) %>% 
  select(c("anno_batch", "sex", "anno_age",  "ADdiag2types", "ADdiag3types" , "brainRegion", "study"))

colnames(NPH_pateint) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
rownames(NPH_pateint) = NULL

saveRDS(NPH_pateint, "/cosmos/data/downloaded-data/AD-meta/Kamath-NPH/data/0-source/meta/Kamath-NPH-patient.rds" )



### Kamath- other ===========================================================================================================================
kamath_dir = "/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/Kamath/annotations/rds-list"
kamath_data <- readRDS(file.path(kamath_dir, "kamath-meta.rds"))

Human_meta <- bind_rows(lapply(kamath_data, function(df) { 
  dplyr::filter(df, anno_organism == "human" & ds_batch != "human_NPH" & anno_condition %in% c("Ctrl" , "AD", "AD_Braak2", "AD_Braak6"))
}))


datasets <- unique(Human_meta$ds_batch)

meta <- list()
for (data in datasets) {
  meta[[data]] = filter(Human_meta, ds_batch == data)
  print(data)
  print(paste("patients", length(unique(meta[[data]]$anno_batch)))) 
  print(paste("cells", length(unique(rownames(meta[[data]]))))) 
  print(paste("dimentions of meta data", dim(meta[[data]])))
  print("next")
  
}


other_pateint <- Human_meta %>% select(c("anno_batch", "anno_sex", "anno_age",  "anno_condition", "anno_region")) %>% unique() %>% 
  mutate(sex = ifelse(anno_sex == "M", "male", "female"), 
         study = "Kamath-other", 
         brainRegion = anno_region,  
         ADdiag2types = ifelse(anno_condition == "Ctrl", "no", "yes"), 
         ADdiag3types = anno_condition) %>% 
  select(c("anno_batch", "sex", "anno_age",  "ADdiag2types", "ADdiag3types" , "brainRegion", "study"))

colnames(NPH_pateint) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")
rownames(NPH_pateint) = NULL

saveRDS(NPH_pateint, "/cosmos/data/downloaded-data/AD-meta/Kamath-NPH/data/0-source/meta/Kamath-NPH-patient.rds" )




### ling - 2024 ===========================================================================================================================
meta_patient <- data.table::fread("/space/scratch/nairuz-rewiring-normal/data/raw/ling/SZvillage_donorMetadata.txt")
meta_patient_final <- meta_patient %>%  mutate(sex = ifelse (Sex == "Female", "female", "male"),
                                               brainRegion = "dlPFC_B46", 
                                               study = "ling-2024", 
                                               ADdiag3types = NA ,
                                               ADdiag2types = ifelse(Schizophrenia == "Unaffected", "no", "yes"))%>% 
  select(c("nbb.id", "sex", "Age", "ADdiag2types", "ADdiag3types" , "brainRegion", "study")) 


colnames(meta_patient_final) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")

saveRDS(meta_patient_final, "/space/scratch/nairuz-rewiring/data/0.processed-raw/all-patient-metaData.rds")


### Ruzicka - 2024 ===========================================================================================================================
ref <- readRDS("/space/scratch/nairuz-rewiring/data/0.processed-raw/all-patient-metaData.rds")

cmc_patient <- readRDS("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/CMC/CMC-patients.rds")
cmc_patient_final <- cmc_patient %>%  mutate(
                                               study = "CMC-MSSM-2024", 
                                               ADdiag3types = NA ,
                                               ADdiag2types = disorder )%>% 
  select(c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types" , "brainRegion", "study")) 



multi_patient <- readRDS("/space/scratch/rui_sz_project/PavlabSZProject/data/data_processed/SZBDMulti-Seq/SZBDMulti-Seq-patients.rds")
multi_patient_final <- multi_patient %>%  mutate(
  study = "SZBDMulti-Seq-2024", 
  ADdiag3types = NA ,
  ADdiag2types = disorder )%>% 
  select(c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types" , "brainRegion", "study")) 


saveRDS(rbind(ref, cmc_patient_final,multi_patient_final), "/space/scratch/nairuz-rewiring/data/0.processed-raw/all-patient-metaData.rds")


### re-update all patients ===========================================================================================================================
all_patients <- readRDS("/space/scratch/nairuz-rewiring/data/1.processed-5per/all-patient-metaData.rds") # Load the existing all_patients file

# Load the new patient metadata files
kamath <- readRDS("/cosmos/data/downloaded-data/AD-meta/Z-META-DATA/total-patients/Kamath-NPH-patient.rds")
deJager <-  readRDS("/cosmos/data/downloaded-data/AD-meta/Z-META-DATA/total-patients/ROSMAP-DeJager-patient.rds")
erosion <- readRDS("/cosmos/data/downloaded-data/AD-meta/Z-META-DATA/total-patients/ROSMAP-erosion-patient.rds")
micro <- readRDS("/cosmos/data/downloaded-data/AD-meta/Z-META-DATA/total-patients/ROSMAP-Micro-patient.rds")
MIT <- readRDS("/cosmos/data/downloaded-data/AD-meta/Z-META-DATA/total-patients/ROSMAP-MIT-patient.rds")

new_patients <- do.call(rbind, list(kamath, deJager, erosion, micro, MIT)) # Combine all new patient data into one dataframe
all_patients_updated <- rbind(all_patients, new_patients) # Combine new patients with the existing all_patients data


# Sanity check
print(paste("New patients dimensions: ", dim(new_patients)))
print(paste("All patients dimensions (before): ", dim(all_patients)))
print(paste("All patients dimensions (after): ", dim(all_patients_updated)))

# Save the updated all_patients file
saveRDS(all_patients_updated, "/space/scratch/nairuz-rewiring/data/1.processed-5per/all-patient-metaData.rds")


### re-update all patients ===========================================================================================================================
all_patients <- readRDS("/space/scratch/nairuz-rewiring/data/1.processed-5per/all-patient-metaData.rds") # Load the existing all_patients file
mathys1 <- readRDS("/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/meta/ROSMAP-Mathys1-patient.rds") 
ROSMAP_clinical <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/MIT-ROSMAP/meta/ROSMAP_clinical.csv")[ , c("projid", "msex", "age_death", "Study")] %>% 
  filter(Study %in% c("ROS", "MAP")) #mathys1? 
ROSMAP_clinical$patientID = ROSMAP_clinical$projid

mathys1_updated <- left_join(mathys1, ROSMAP_clinical[, c("patientID", "age_death")])
mathys1_updated <- mathys1_updated %>%  mutate(age = age_death) %>% 
  select(c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types" , "brainRegion", "study")) 

all_patients_updated <- rbind(all_patients, mathys1_updated) # Combine new patients with the existing all_patients data








### figuring out patient overlap
individual_metadata <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/MIT-ROSMAP/meta/individual_metadata.csv")[, c("individualID", "sex", "subject")]
multiome_metadata <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/MIT-ROSMAP/meta/multiome_metadata.csv")[, c("specimenID", "assay")]
biospecimen_metadata <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/MIT-ROSMAP/meta/biospecimen_metadata.csv")[, c("individualID", "specimenID", "tissue")] 







### SEA-AD- MTG ===========================================================================================================================
SEA_AD_MTG <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/meta/SEA-AD-MIT-cells-meta.csv")
colnames(SEA_AD_MTG) <- gsub(" ", "_", colnames(SEA_AD_MTG))

SEAAD_MTG_patient <- SEA_AD_MTG[, c("Donor_ID", "Sex", "Age_at_Death", "Overall_AD_neuropathological_Change", "Cognitive_Status", "Brain_Region")] %>% mutate(study = "SEA-AD") 
SEAAD_MTG_patient <- unique(SEAAD_MTG_patient)

SEAAD_MTG_patient <- SEAAD_MTG_patient[ !(Cognitive_Status == "Reference"), ]


SEAAD_MTG_patient <- SEAAD_MTG_patient %>% 
  mutate(sex = ifelse (Sex == "Female", "female", "male"), 
         study = "SEA-AD-MTG",
         brainRegion = Brain_Region, 
         ADdiag2types = "yes", 
         ADdiag3types =  Overall_AD_neuropathological_Change ) %>% 
  select(c("Donor_ID", "sex", "Age_at_Death", "ADdiag2types", "ADdiag3types" , "Cognitive_Status", "brainRegion", "study"))

colnames(SEAAD_MTG_patient) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "CognitiveStatus",  "brainRegion", "study")
saveRDS(SEAAD_MTG_patient, "/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/meta/SEA-AD-MTG-patient.rds" )



### SEA-AD- DLPFC   ===========================================================================================================================
SEA_AD_DLPFC <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/SEA-AD/DLPFC/RNAseq/SEAAD_DLPFC_RNAseq_final_meta.csv")
colnames(SEA_AD_DLPFC) <- gsub(" ", "_", colnames(SEA_AD_DLPFC))

SEAAD_DLPFC_patient <- SEA_AD_DLPFC[, c("Donor_ID", "Sex", "Age_at_Death", "Overall_AD_neuropathological_Change", "Cognitive_Status", "Brain_Region")] %>% mutate(study = "SEA-AD") 
SEAAD_DLPFC_patient <- unique(SEAAD_DLPFC_patient)

SEAAD_DLPFC_patient <- SEAAD_DLPFC_patient[ !(Cognitive_Status == "Reference"), ]


SEAAD_DLPFC_patient <- SEAAD_DLPFC_patient %>% 
  mutate(sex = ifelse (Sex == "Female", "female", "male"), 
         study = "SEA-AD-MTG",
         brainRegion = Brain_Region, 
         ADdiag2types = "yes", 
         ADdiag3types =  Overall_AD_neuropathological_Change ) %>% 
  select(c("Donor_ID", "sex", "Age_at_Death", "ADdiag2types", "ADdiag3types" , "Cognitive_Status", "brainRegion", "study"))

colnames(SEAAD_DLPFC_patient) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "CognitiveStatus",  "brainRegion", "study")
saveRDS(SEAAD_DLPFC_patient, "/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/meta/SEA-AD-DLPFC-patient.rds" )









SEA_DLPFC <- data.table::fread("/cosmos/data/downloaded-data/AD-nairuz/test/CELL-issue/SEA-AD/DLPFC/RNAseq/SEAAD_DLPFC_RNAseq_final_meta.csv")
colnames(SEA_DLPFC) <- gsub(" ", "_", colnames(SEA_DLPFC))

SEA_patient_DLPFC <- SEA_DLPFC[, c("Donor_ID", "Sex", "Age_at_Death", "Overall_AD_neuropathological_Change", "Brain_Region")] %>% mutate(study = "SEA-AD") 
SEA_patient_DLPFC <- unique(SEA_patient_DLPFC)
colnames(SEA_patient_DLPFC) <- c("patientID", "sex", "age", "diagnosis", "brainRegion", "study")














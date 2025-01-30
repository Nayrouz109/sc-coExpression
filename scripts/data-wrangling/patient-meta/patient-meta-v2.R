
#### Nov 13th ===============================================================================================================
#### Nov 13th ===============================================================================================================
#### Nov 13th ===============================================================================================================
all_patients <- read_csv("/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv")

all_patients <- all_patients %>%

  rename(diagnosis = ADdiag2types) %>% 
  mutate(study = case_when(
    study == "ROSMAP-DeJager" ~ "ROSMAP-DeJager-2024",
    study == "ROSMAP-MIT" ~ "ROSMAP-Mathys-2023",
    TRUE ~ study
  )) %>%
  
  mutate(disease = if_else(study %in% c("ling-2024", "CMC-MSSM-2024", "SZBDMulti-Seq-2024"), 
                           "SCZ", "AD")) %>%
  
  rename(patientID_brainRegion_disease = patientID_brainRegion) %>% 
  mutate(patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis))



#### psychAD ===============================================================================================================
psychAD_patients <- read_csv("/cosmos/data/downloaded-data/AD-meta/psychAD-2024/data/meta/patient-metaData-final.csv")

MSBB <- psychAD_patients %>% filter(cohort == "MSBB")
HBCC <- psychAD_patients %>% filter(cohort == "HBCC")
RADC <- psychAD_patients %>% filter(cohort == "RADC")

# Transform MSBB dataset as given
MSBB <- MSBB %>%  
  mutate(study = "MSBB-2024", 
         brainRegion = "dlPFC", 
         disease = diagnosis, 
         diagnosis = ifelse(disease == "control", "no", "yes"), 
         patientID_brainRegion_disease = paste0(DonorID, "_", brainRegion, "_", diagnosis), 
         ADdiag3types = diagnosis) %>% 
  select(c("DonorID", "sex", "ageDeath", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease"))

# Transform HBCC dataset
HBCC <- HBCC %>%  
  mutate(study = "HBCC-2024", 
         brainRegion = "dlPFC", 
         disease = diagnosis, 
         diagnosis = ifelse(disease == "control", "no", "yes"), 
         patientID_brainRegion_disease = paste0(DonorID, "_", brainRegion, "_", diagnosis), 
         ADdiag3types = diagnosis) %>% 
  select(c("DonorID", "sex", "ageDeath", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease"))

# Transform RADC dataset
RADC <- RADC %>%  
  mutate(study = "RADC-2024", 
         brainRegion = "dlPFC", 
         disease = diagnosis, 
         diagnosis = ifelse(disease == "control", "no", "yes"), 
         patientID_brainRegion_disease = paste0(DonorID, "_", brainRegion, "_", diagnosis), 
         ADdiag3types = diagnosis) %>% 
  select(c("DonorID", "sex", "ageDeath", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease"))

# Rename columns in all three datasets
colnames(MSBB) <- c("patientID", "sex", "age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")
colnames(HBCC) <- c("patientID", "sex", "age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")
colnames(RADC) <- c("patientID", "sex", "age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")

all_patientsV1 <- bind_rows(all_patients, MSBB, HBCC, RADC)


#### hoffman ===============================================================================================================

hoffman <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/Hoffman-2023/data/expr/patient_metadata") %>% unique()

hoffman <- hoffman %>%  
  mutate(study = "Hoffman-2023", 
         sex = ifelse (Sex == "Female", "female", "male"), 
         brainRegion = "dlPFC", 
         disease = "AD", 
         diagnosis = ifelse(Dx_AD == "Control", "no", "yes"), 
         patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis), 
         ADdiag3types = Dx_AD) %>% 
  select(c("patientID", "sex", "Age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")) %>% unique()

# Rename columns in all three datasets
colnames(hoffman) <- c("patientID", "sex", "age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")

all_patientsV1 <- all_patientsV1 %>% select(-1)
all_patientsV2 <- rbind(all_patientsV1, hoffman)



#### TREM ===============================================================================================================
TREM <- readRDS("/cosmos/data/downloaded-data/AD-meta/ROSMAP-TREM2/data/1-filt/ROSMAP-TREM2-2020-patient.rds")
TREM <- TREM %>% 
  rename(diagnosis = ADdiag2types) %>% 
  mutate(disease = "AD") %>% 
  mutate(patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis))

all_patientsV3 <- rbind(all_patientsV2, TREM)


#### mathys 2019 ===============================================================================================================
dir = "/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/meta"
snRNAseqPFC_BA10_biospecimen_metadata    <- read_csv(file.path(dir, "snRNAseqPFC_BA10_biospecimen_metadata.csv")) %>% 
  rename(patientID = projid)

snRNAseqPFC_BA10_biospecimen_metadata$patientID = as.character(snRNAseqPFC_BA10_biospecimen_metadata$patientID)
ROSMAP_clinical_mathys2019 <- read_csv(file.path(dir, "ROSMAP_clinical.csv")) %>% 
  rename(patientID = projid)

mathys2019p.V0 <- left_join(snRNAseqPFC_BA10_biospecimen_metadata, ROSMAP_clinical_mathys2019, by = "patientID")   
mathys2019p.V0 <- mathys2019p.V0 %>% select(c("patientID", "individualID.x", "tissue", "BrodmannArea", "msex", 
                                              "apoe_genotype", "age_death", "pmi", "braaksc", "ceradsc", "cogdx", "dcfdx_lv")) 



mathys2019p = readRDS("/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/meta/ROSMAP-Mathys1-patient.rds")
mathys2019p$patientID = as.character(mathys2019p$patientID)
mathys2019p <- left_join(mathys2019p, mathys2019p.V0[, c("patientID", "age_death")] , by = "patientID" )

mathys2019p <- mathys2019p %>% 
  rename(age = age_death) %>% 
  rename(diagnosis = ADdiag2types) %>% 
  mutate(disease = "AD") %>% 
  mutate(patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis)) %>% 
  mutate(brainRegion = "dlPFC_BA10") %>% 
  mutate(study = "ROSMAP-Mathys-2019")

saveRDS(mathys2019p.V0, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-patient-metaData/snRNAseqPFC_BA10/mathys2019-wrangled.rds")
all_patientsV4 <- rbind(all_patientsV3, mathys2019p)



#### mathys 2024 ===============================================================================================================

multi_region <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT-multiregion-2024/data/all_brain_regions_filt_preprocessed_scanpy_norm.final_noMB.cell_labels.tsv") %>% select(c("region", "projid")) %>% unique()
ROSMAP_clinical <- read_csv("/cosmos/data/downloaded-data/AD-meta/ROSMAP-Mathys1/data/0-source/meta/ROSMAP_clinical.csv") %>% select(c("projid", "apoe_genotype", "age_death", 
                                                                                                                                       "braaksc", "ceradsc", "cogdx", "dcfdx_lv", "individualID"))
multi_region$projid    = as.character(multi_region$projid) 
ROSMAP_clinical$projid = as.character(ROSMAP_clinical$projid) 
multi_region.V0 <- left_join(multi_region, ROSMAP_clinical, by = "projid") 
multi_region.V0 <- multi_region.V0 %>% rename(patientID = projid) %>% select(-("region"))

MIT_dir = "/cosmos/data/downloaded-data/AD-meta/ROSMAP-patient-metaData/MIT-multiome"
MIT_patients <- read_csv(file.path(MIT_dir, "MIT_ROSMAP_Multiomics_individual_metadata.csv")) %>% select(c("individualID", "sex", "subject")) %>% filter(subject %in% diagnosis$subject) %>% unique()
diagnosis <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT-multiregion-2024/data/MIT-multiome") %>% select(-c(age_death))
diagnosis <- diagnosis %>%
  left_join(MIT_patients, by = "subject")

multi_region.V1 <- multi_region.V0 %>% left_join(diagnosis, by = "individualID") %>% unique()


# Update specific region values with corresponding BA regions
multi_region.V1 <- multi_region.V1 %>%
  mutate(region = case_when(
    region == "AG"  ~ "AG_BA39.40",
    region == "MT"  ~ "MT_BA22",
    region == "PFC" ~ "PFC_BA10",
    region == "EC"  ~ "EC_BA28",
    TRUE            ~ region  # Keep the existing value if it doesn't match any specified case
  ))


saveRDS(multi_region.V1, "/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT-multiregion-2024/data/annData-cellType/mathys2024-wrangled.rds")


multi_region.V2 <- multi_region.V1 %>% 
  rename(age = age_death) %>% 
  mutate(ADdiag3types = clinicalDiag) %>% 
  mutate(diagnosis = ifelse(pathAD == "non-AD", "no", "yes")) %>% 
  mutate(disease = "AD") %>%
  mutate(brainRegion = region) %>% 
  mutate(patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis)) %>% 
  mutate(study = "ROSMAP-Mathys-2024") %>% 
  select(c("patientID", "sex", "age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")) %>% unique()


colnames(hoffman) <- c("patientID", "sex", "age", "diagnosis", "ADdiag3types", "brainRegion", "study", "patientID_brainRegion_disease", "disease")
all_patientsV5 <- rbind(all_patientsV4, multi_region.V2)


##### seaAD =========================================================================================================
seaAD_MTG   <- readRDS("/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/expr/annData-cellTpe/00sea_MTG_patient.rds")
seaAD_DLPFC <- readRDS("/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/annData-cellTpe/00sea_DLPFC_patient.rds")

seaAD_MTG.V1 <- seaAD_MTG %>% 
  rename(patientID = Donor_ID, age = Age_at_Death) %>% 
  mutate(
    sex = ifelse(Sex == "Female", "female", "male"),
    diagnosis = ifelse(nairuz_category == "control", "no", "yes"),
    ADdiag3types = Overall_AD_neuropathological_Change,
    brainRegion = Brain_Region,
    study = "SEA-AD-MTG-2024",
    patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis),
    disease = "AD"
  )


seaAD_DLPFC.V1 <- seaAD_DLPFC %>% 
  rename(patientID = Donor_ID, age = Age_at_Death) %>% 
  mutate(
    sex = ifelse(Sex == "Female", "female", "male"),
    diagnosis = ifelse(nairuz_category == "control", "no", "yes"),
    ADdiag3types = Overall_AD_neuropathological_Change,
    brainRegion = Brain_Region,
    study = "SEA-AD-DLPFC-2024",
    patientID_brainRegion_disease = paste0(patientID, "_", brainRegion, "_", diagnosis),
    disease = "AD"
  )




combined_seaAD.v1 <- combined_seaAD %>% select(colnames(all_patientsV5)) 
all_patientsV6 <- rbind(all_patientsV5, combined_seaAD.v1)

### sanity check =====================================
# make sure study names are consistent between file names and study names in patient metadata file 
matrices_dir = "/space/scratch/nairuz-rewiring/data/0.prcsd-raw"
matrices <- list.files(matrices_dir, pattern = "h5ad")
study_names <-unique(sub("^[^_]+_(.*)\\.h5ad$", "\\1", matrices))

#> study_names[(!(study_names %in% all_patientsV6$study))]
#[1] "Kamath-NPH-2023"     "MSSM-2024"           "ROSMAP-erosion-2023" "SZBDMulti-2024"    

all_patientsV6 <- all_patientsV6 %>%
  mutate(study = case_when(
    study == "Kamath-NPH" ~ "Kamath-NPH-2023",
    study == "ROSMAP-erosion" ~ "ROSMAP-erosion-2023",
    study == "SZBDMulti-Seq-2024" ~ "SZBDMulti-2024", 
    study == "CMC-MSSM-2024" ~ "MSSM-2024" , 
    TRUE ~ study
  ))

length(unique(all_patientsV6$study))
length(unique(study_names))

all(study_names %in% all_patientsV6$study)

write.csv(all_patientsV6, "/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv", row.names = TRUE)
studies_p = ((data.table::fread("/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv", header = TRUE))$study)
all(study_names %in% studies_p)

length((filter(all_patientsV6, diagnosis == "no"))$patientID_brainRegion_disease)
length((filter(all_patientsV6, diagnosis == "yes"))$patientID_brainRegion_disease)








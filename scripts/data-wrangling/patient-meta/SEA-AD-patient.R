

##### MTG ====================================================================================================================
sea_MTG_cells <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/meta/SEA-AD-MIT-cells-meta.csv") 
#sea_MTG       <- readRDS("/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/meta/SEA-AD-MTG-patient.rds")


metadata_columns <- c(
  # Demographic Information
  "Donor ID", "Sex", "Gender", "Age at Death", 
  #"Race (choice=White)", "Race (choice=Black/ African American)", 
  #"Race (choice=Asian)", "Race (choice=American Indian/ Alaska Native)", 
  #"Race (choice=Native Hawaiian or Pacific Islander)", "Race (choice=Unknown or unreported)", 
  #"Race (choice=Other)", "Hispanic/Latino",
  
  # Clinical and Health Information
  "PMI", "Overall AD neuropathological Change", "Thal", "Braak", 
  "CERAD score", "Overall CAA Score", "Highest Lewy Body Disease", 
  "Total Microinfarcts (not observed grossly)", "Total microinfarcts in screening sections", 
  "Atherosclerosis", "Arteriolosclerosis", "LATE", "Neurotypical reference" , 
  
  # Cognitive and Diagnostic Information
  "Cognitive Status", #"Last CASI Score", "Interval from last CASI in months", 
  #"Last MMSE Score", "Interval from last MMSE in months", 
  #"Last MOCA Score", "Interval from last MOCA in months", 
  "APOE Genotype",
  
  # Study and Sample Information
  "Primary Study Name", "Secondary Study Name", "Brain Region"
)

patient_metadata <- sea_MTG_cells[, ..metadata_columns] %>% unique()

# Control Group: Neurotypical or no AD pathology
control_cells <- patient_metadata %>%
  filter(`Neurotypical reference` == TRUE | 
           `Overall AD neuropathological Change` == "Not AD") %>% 
  mutate(nairuz_category = "control")

# Alzheimer's Group: AD pathology based on neuropathological assessment
alzheimers_cells <- patient_metadata %>%
  filter(`Overall AD neuropathological Change` %in% c("Intermediate", "High", "Low")) %>% 
  mutate(nairuz_category = "earlyAD")

# Advanced AD cases
advanced_ad_cells <- alzheimers_cells %>%
  filter(Braak %in% c("Braak V", "Braak VI"),
         Thal %in% c("Thal 4", "Thal 5"),
         `CERAD score` %in% c("Moderate", "Frequent")) %>% 
  mutate(nairuz_category = "lateAD")

sea_MTG_patient <- rbind(control_cells, advanced_ad_cells)
sea_MTG_patient <- rbind(sea_MTG_patient, alzheimers_cells %>% filter (!(`Donor ID` %in% sea_MTG_patient$`Donor ID` )))

colnames(sea_MTG_patient) <- gsub(" ", "_", colnames(sea_MTG_patient))
sea_MTG_patient <- (sea_MTG_patient[,  lapply(.SD, function(x) gsub(" ", "_", as.character(x)))])


saveRDS(sea_MTG_patient, "/cosmos/data/downloaded-data/AD-meta/SEA-AD-MTG/0-source/expr/annData-cellTpe/00sea_MTG_patient.rds")

##### dlPFC ====================================================================================================================
#SEAAD_DLPFC_patient <- readRDS("/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/meta/SEA-AD-DLPFC-patient.rds") 
sea_DLPFC_cells <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/SEAAD_A9_RNAseq_final-nuclei_metadata.2024-02-13.csv")

colnames(sea_DLPFC_cells) <- gsub(" ", "_", colnames(sea_DLPFC_cells))
metadata_columns <- c(
  # Demographic Information
  "Donor_ID", "Sex", "Gender", "Age_at_Death", 
  # "Race_(choice=White)", "Race_(choice=Black/_African_American)", 
  # "Race_(choice=Asian)", "Race_(choice=American_Indian/_Alaska_Native)", 
  # "Race_(choice=Native_Hawaiian_or_Pacific_Islander)", "Race_(choice=Unknown_or_unreported)", 
  # "Race_(choice=Other)", "Hispanic/Latino",
  
  # Clinical and Health Information
  "PMI", "Overall_AD_neuropathological_Change", "Thal", "Braak", 
  "CERAD_score", "Overall_CAA_Score", "Highest_Lewy_Body_Disease", 
  "Total_Microinfarcts_(not_observed_grossly)", "Total_microinfarcts_in_screening_sections", 
  "Atherosclerosis", "Arteriolosclerosis", "LATE", "Neurotypical_reference", 
  
  # Cognitive and Diagnostic Information
  "Cognitive_Status", # "Last_CASI_Score", "Interval_from_last_CASI_in_months", 
  # "Last_MMSE_Score", "Interval_from_last_MMSE_in_months", 
  # "Last_MOCA_Score", "Interval_from_last_MOCA_in_months", 
  "APOE_Genotype",
  
  # Study and Sample Information
  "Primary_Study_Name", "Secondary_Study_Name", "Brain_Region"
)


patient_metadata <- sea_DLPFC_cells[, ..metadata_columns, drop = FALSE] %>% unique()
patient_metadata <- patient_metadata[, lapply(.SD, function(x) gsub(" ", "_", as.character(x))), .SDcols = metadata_columns] # Replace spaces in values within the selected columns

# Define the control group based on Neurotypical reference or no AD pathology
control_cells <- patient_metadata %>%
  filter(Neurotypical_reference == TRUE | Overall_AD_neuropathological_Change == "Not_AD") %>% 
  mutate(nairuz_category = "control")

# Define the Alzheimer's group based on neuropathological assessment
alzheimers_cells <- patient_metadata %>%
  filter(Overall_AD_neuropathological_Change %in% c("Intermediate", "High", "Low")) %>% 
  mutate(nairuz_category = "earlyAD")

# Define the advanced AD cases group based on specific criteria
advanced_ad_cells <- alzheimers_cells %>%
  filter(Braak %in% c("Braak_V", "Braak_VI"),
         Thal %in% c("Thal_4", "Thal_5"),
         CERAD_score %in% c("Moderate", "Frequent")) %>% 
  mutate(nairuz_category = "lateAD")

# Combine groups, ensuring no duplication of Donor IDs
sea_DLPFC_patient <- rbind(control_cells, advanced_ad_cells)
sea_DLPFC_patient <- rbind(sea_DLPFC_patient, alzheimers_cells %>% filter(!(`Donor_ID` %in% sea_DLPFC_patient$`Donor_ID`)))

length(unique((filter(sea_DLPFC_patient, nairuz_category == "control"))$Donor_ID))
length(unique((filter(sea_DLPFC_patient, Overall_AD_neuropathological_Change == "Not_AD"))$Donor_ID))
length(unique((filter(sea_DLPFC_patient, Overall_AD_neuropathological_Change == "Reference"))$Donor_ID))
length(unique((filter(sea_DLPFC_patient, nairuz_category == "lateAD"))$Donor_ID))
length(unique((filter(sea_DLPFC_patient, nairuz_category == "earlyAD"))$Donor_ID))


saveRDS(sea_DLPFC_patient, "/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/annData-cellTpe/00sea_DLPFC_patient.rds")



### how i decided on categorizing AD stages =======================================================================================
length(unique(sea_MTG_cells$`Neurotypical reference`))   ***
FALSE  TRUE

(unique(sea_MTG_cells$`Overall AD neuropathological Change`))   ***
"Intermediate" "High"         "Low"          "Not AD"       "Reference"   



> (unique(sea_MTG_cells$Thal))
[1] "Thal 3"    "Thal 5"    "Thal 2"    "Thal 4"    "Thal 0"    "Reference" "Thal 1"  

> (unique(sea_MTG_cells$Braak))
[1] "Braak IV"  "Braak V"   "Braak VI"  "Braak II"  "Braak III" "Reference" "Braak 0"  

> (unique(sea_MTG_cells$`CERAD score`))
[1] "Absent"    "Moderate"  "Sparse"    "Frequent"  "Reference"

> (unique(sea_MTG_cells$`Overall CAA Score`))
[1] "Moderate"       "Mild"           "Not identified" "Severe"         "Reference"  

> (unique(sea_MTG_cells$`Highest Lewy Body Disease`))
[1] "Not Identified (olfactory bulb not assessed)" "Not Identified (olfactory bulb assessed)"     "Olfactory bulb only"                          "Limbic (Transitional)"                       
[5] "Neocortical (Diffuse)"                        "Amygdala-predominant"                         "Reference"                                    "Brainstem-predominant"   

> (unique(sea_MTG_cells$LATE))
[1] "Not Identified" "LATE Stage 2"   "LATE Stage 1"   "LATE Stage 3"   "Reference"      "Unclassifiable"

> (unique(sea_MTG_cells$`Cognitive Status`))
[1] "No dementia" "Dementia"    "Reference"

> (unique(sea_MTG_cells$`APOE Genotype`))
[1] "2/2" "3/3" "3/4" "2/3" "4/4" ""    "2/4"

> length(unique(sea_MTG_cells$`Donor ID`))
[1] 89

> (unique(sea_MTG_cells$`Brain Region`))
[1] "Human MTG"            "Human MTG_L5"         "Human MTG All Layers"






















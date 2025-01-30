

path = "/cosmos/data/downloaded-data/AD-meta/psychAD-2024/data/meta"
media = read_csv(file.path(path, "media-1.csv"))
IDs_mapping = read_csv(file.path(path, "NPS-AD_DonorID_mapping.csv")) 
IDs_mapping$IndividualID = IDs_mapping$individualID
diagnosis_metadata   = read_csv(file.path(path, "NPS-AD_neuropathology_and_diagnosis_metadata.csv"))

patient_metaData <- left_join(media, IDs_mapping)
patient_metaData <- left_join(patient_metaData, diagnosis_metadata, by = "IndividualID" )



patient_metaDataV <- patient_metaData %>%
  mutate(diagnosis = case_when(
    crossDis_AD == 1 ~ "AD",
    dx_AD == 1 ~ "AD",
    crossDis_SCZ == 1 ~ "SCZ",
    crossDis_DLBD == 1 ~ "DLBD",
    crossDis_Vas == 1 ~ "Vas",
    crossDis_BD == 1 ~ "BD",
    crossDis_Tau == 1 ~ "Tau",
    crossDis_PD == 1 ~ "PD",
    crossDis_FTD == 1 ~ "FTD",
    PD == 0 & SCZ == 0 & BD == 0 & 
      (is.na(FTD) | FTD == 0) &
      (is.na(MS) | MS == 0) &
      (is.na(PTSD) | PTSD == 0) &
      (is.na(ADHD) | ADHD == 0) &
      (is.na(OCD) | OCD == 0) &
      (is.na(dx_AD) | dx_AD == 0) &
      (is.na(crossDis_DLBD) | crossDis_DLBD == 0) &
      (is.na(crossDis_Vas) | crossDis_Vas == 0) &
      (is.na(crossDis_Tau) | crossDis_Tau == 0) &
      (is.na(CERAD.x) | CERAD.x == 1) &
      (is.na(CDR) | CDR == 0) &
      (is.na(BRAAK_AD) | BRAAK_AD %in% c(0, 1, 2, 3)) ~ "control",
    TRUE ~ "other"
  ))

summary <- patient_metaDataV %>%
  group_by(cohort, diagnosis) %>%
  summarize(total_patients = n(), .groups = "drop")


write.csv(patient_metaDataV, "/cosmos/data/downloaded-data/AD-meta/psychAD-2024/data/meta/patient-metaData-final.csv")




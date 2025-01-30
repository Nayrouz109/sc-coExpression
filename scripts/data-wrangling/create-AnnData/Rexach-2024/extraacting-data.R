

### read in the data ================================
#meta <- data.table::fread("/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/liger_subcluster_metadata_v2.csv")
rds <- readRDS("/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/pci_seurat.rds")


### save patient meta data ===================================================
all_patient_metaData <- rds@meta.data %>%
  select(library_id, autopsy_id, finalsite, clinical_dx, npdx1, sex, region, library_batch, seq_batch,sample, 
         age, pmi, ) %>%
  filter(clinical_dx %in% c("Control", "AD")) %>% unique()
rownames(all_patient_metaData) <- NULL

all_patients <- readRDS("/space/scratch/nairuz-rewiring/data/1.processed-5per/all-patient-metaData.rds")

rexach_p <- all_patient_metaData %>% 
  mutate(sex = ifelse (sex == "F", "female", "male"), 
         study = "Rexach-2024",
         brainRegion = region, 
         ADdiag2types = ifelse(clinical_dx == "AD", "yes", "no"), 
         ADdiag3types =  ADdiag2types ) %>% 
  select(c("autopsy_id", "sex", "age", "ADdiag2types", "ADdiag3types" , "brainRegion", "study"))

rownames(rexach_p) <- NULL
colnames(rexach_p) = c("patientID", "sex", "age", "ADdiag2types", "ADdiag3types", "brainRegion", "study")

length(filter(rexach_p, ADdiag2types == "yes")$patientID) #27
length(filter(rexach_p, ADdiag2types == "no")$patientID)  #25 

all_patients_updated <- rbind(all_patients, rexach_p)
all_patients_updated$patientID_brainRegion <- paste0(all_patients_updated$patientID, all_patients_updated$brainRegion)

write.csv(all_patients_updated, "/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv")
write.csv(all_patient_metaData, "/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/all_patient_metaData.csv")


### Extract Components for annData ===================================================

metaExtracted <- rds@meta.data %>% mutate (patientID = autopsy_id) %>%  filter(clinical_dx %in% c("Control", "AD"))
mtxExtracted  <- rds@assays$RNA@counts[, colnames(rds@assays$RNA@counts) %in% rownames(metaExtracted)]

all_patients_updated <- data.table::fread("/space/scratch/nairuz-rewiring/data/0.prcsd-raw/all-patient-metaData.csv", header = TRUE) 
all(metaExtracted$patientID %in% all_patients_updated$patientID)

Matrix::writeMM(mtxExtracted, "/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components/countsV2.mtx")
Matrix::writeMM(mtxExtracted, gzfile("/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components/countsV2.mtx.gz"))
write.csv(rownames(mtxExtracted), "/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components/genes.csv", row.names = FALSE)
write.csv(colnames(mtxExtracted), "/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components/cell_ids.csv", row.names = FALSE)
write.csv(metaExtracted, "/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components/meta_data.csv", row.names = TRUE)









### exploring the metadata ===================================================
matrixxx <- Matrix::readMM("/cosmos/data/downloaded-data/AD-meta/Rexach-2024/data/0-source/seurat-components/counts.mtx")



class(rds@assays$RNA@counts)
class(rds@assays$RNA@data)

dim(rds@assays$RNA@counts)
dim(rds@assays$RNA@data)

colnames(rds@assays$RNA@counts)[1:3]
colnames(rds@assays$RNA@data)[1:3]
rownames(rds@assays$RNA@counts)[1:3]
rownames(rds@assays$RNA@data)[1:3]

dim(rds@meta.data)


### filter only AD and control samples 
unique(rds@meta.data$autopsy_id)
unique(rds@meta.data$sample)
unique(rds@meta.data$finalsite)
unique(rds@meta.data$clinical_dx) #this one 
unique(rds@meta.data$npdx1) #this one 
unique(rds@meta.data$region)
unique(rds@meta.data$cell_type)
unique(rds@meta.data$cluster_cell_type)


meta_patients <- rds@meta.data[, c("autopsy_id", "sample", "clinical_dx", "npdx1", "region")] %>% unique()

length(unique((meta_patients %>% filter(meta_patients$npdx1 == "PSP"))$autopsy_id))
length(unique((meta_patients %>% filter(meta_patients$npdx1 == "AD"))$autopsy_id))
length(unique((meta_patients %>% filter(meta_patients$npdx1 == "PiD"))$autopsy_id))
length(unique((meta_patients %>% filter(meta_patients$npdx1 == "Normal"))$autopsy_id))
length(unique((meta_patients %>% filter(meta_patients$npdx1 == "Control"))$autopsy_id))
length(unique((meta_patients %>% filter(meta_patients$npdx1 == "Pick's"))$autopsy_id))









PB = readRDS("/cosmos/data/project-data/NW-rewiring/data/1.prcsd-clean/pseudobulk/cpm-log/Mic_MSBB-2024.rds")
PB$meta <- PB$meta %>%  filter(disease %in% c("AD", "CTL"))

PB$expr <- PB$expr[ , colnames(PB$expr) %in% rownames(PB$meta) ]

AD = PB$expr[ , colnames(PB$expr) %in% rownames(PB$meta %>% filter(disease == "AD")) ]
CTL = PB$expr[ , colnames(PB$expr) %in% rownames(PB$meta %>% filter(disease == "CTL")) ]

boxplot(as.matrix(CTL))
boxplot(as.matrix(AD))

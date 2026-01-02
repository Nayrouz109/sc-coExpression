##### read in data ========================================================================================================================
##### =====================================================================================================================================


#### 1. TFs 
TFs = readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS") 
alex_mic_allRank  = readRDS("/home/nelazzabi/rewiring/data/TFs/alex_mic_allRank.rds")
alex_mic_fishZ    = readRDS("/home/nelazzabi/rewiring/data/TFs/alex_mic_fishZ.rds")



ADgene_list <- c(
  "ABCA7", "ABI3", "ACE2", "ADAM10", "APH1B", "APOE", "BIN1", "CD2AP", 
  "CR1", "ECHDC3", "FERMT2", "HS3ST1", "IL34", "MINK1", "PLCG2", "PLD3", 
  "PRKD3", "PTK2B", "SCIMP", "SHARPIN", "SLC2A4", "SORL1", "TREM2", "UNC5C"
)

# Create a list to summarize genes from the SEA-AD MTG paper Feb 11th, 2025 
gene_list <- list(
  inflammatory_processes = c("IL1B", "CSF1R", "STAB1", "NINJ1", "JAK3"),
  interferon_response = c("IRF1", "IRF7", "IFI16"),
  fc_receptors = c("FCGR1A", "FCGR1B", "FCGR2A", "FCGR3B"),
  mhc_class_II_components = c("CD74", "HLA-DRB5"),
  complement_components = c("C1QA", "C1QB"),
  plaque_induced_early = c("CSF1R", "CTSC", "C1QA", "C1QB", "LY86", "FCGR3A"),
  plaque_induced_late = c("CTSD", "CTSS", "LYZ", "APOE"),
  transcription_factors = c("RUNX1", "IKZF1", "NFATC2", "MAF")
)

SEAAD_df <- do.call(rbind, lapply(names(gene_list), function(category) {
  data.frame(Gene = gene_list[[category]], Category = category, stringsAsFactors = FALSE)
}))

eRegulon = fread("/home/nelazzabi/rewiring/data/TFs/eRegulon_metadata_filtered.csv")
eRegulonV1 = eRegulon %>% filter(TF %in% gene_list$transcription_factors) %>% 
  dplyr::select(c("TF", "Gene")) %>%  unique() 

freq = as.data.frame(table(eRegulonV1$Gene))
length(unique((freq %>% filter(Freq > 1))$Var1))






#### 2. my data 
my_data_TF_top200 = readRDS("/cosmos/data/project-data/NW-rewiring/coExpr/disease/TF-targets/tf_top200_list.rds")
my_data_TF_allTrg = readRDS("/cosmos/data/project-data/NW-rewiring/coExpr/disease/TF-targets/tf_allTrg_list.rds") 



my_data_TF_top200_AD  <- my_data_TF_top200[grepl("_AD", names(my_data_TF_top200))]
my_data_TF_allTrg_AD  <- my_data_TF_allTrg[grepl("_AD", names(my_data_TF_allTrg))]
my_data_TF_top200_CTL <- my_data_TF_top200[grepl("_CTL", names(my_data_TF_top200))]
my_data_TF_allTrg_CTL <- my_data_TF_allTrg[grepl("_CTL", names(my_data_TF_allTrg))]



### sanity checks ========================================================================
length(my_data_TF_allTrg_CTL) # should be total # of studies 
length(my_data_TF_top200_CTL) # should be total # of studies 

for (i in 1:length(my_data_TF_allTrg_CTL)) {
  print(length(unique(my_data_TF_allTrg_CTL[[i]]$geneA))) #should all be the same 861TFs 
  #print(length(unique(filtered_data_CTL[[i]]$geneB))) #this is not really that important 
  print(dim(my_data_TF_allTrg_CTL[[i]])) #should all have the same dim (861 * 10,907)
}

for (i in 1:length(my_data_TF_top200_CTL)) {
  print(length(unique(my_data_TF_top200_CTL[[i]]$geneA))) #should all be the same 861TFs 
  #print(length(unique(filtered_data_CTL[[i]]$geneB))) #this is not really that important 
  print(dim(my_data_TF_top200_CTL[[i]])) #should all have the same dim (861 * 200)
}

View(table(my_data_TF_top200_CTL$`Mic_Hoffman-2023_CTL`$geneA))
View(table(my_data_TF_allTrg_CTL$`Mic_Hoffman-2023_CTL`$geneA))

View(table(test_TF_allTargets$geneA)) #optional 
View(table(test_TF_topTargets$geneA)) #optional 

length(unique(alex_mic_allRank$geneA))
length(unique(alex_mic_fishZ$geneA))


TFmydata_vs_micAllRank = unique(my_data_TF_allTrg_CTL$`Mic_Hoffman-2023_CTL`$geneA)[unique(my_data_TF_allTrg_CTL$`Mic_Hoffman-2023_CTL`$geneA) %in% unique(alex_mic_allRank$geneA)]
TFmydata_vs_micFishZ   = unique(my_data_TF_allTrg_CTL$`Mic_Hoffman-2023_CTL`$geneA)[unique(my_data_TF_allTrg_CTL$`Mic_Hoffman-2023_CTL`$geneA) %in% unique(alex_mic_fishZ$geneA)]





#  1. my analysis =======================================================================================================================
combined_df_top200 <- bind_rows(my_data_TF_top200_CTL, .id = "study")
combined_df_alltrg <- bind_rows(my_data_TF_allTrg_CTL, .id = "study")



## 1.1 bar plots - global TF-200 pairs view ***********************************************
tf_target_counts <- combined_df_top200 %>%
  group_by(geneA, geneB) %>%
  summarize(study_count = n(), .groups = "drop") 

### bar plot 
#ggplot(tf_target_counts %>% filter(study_count >= 5), aes(x = study_count)) +
ggplot(tf_target_counts, aes(x = study_count)) +
  geom_bar(fill = "blue", alpha = 0.8) +
  labs(title = "Distribution of TF-top 200 co-expr Pairs Across Studies",
       x = "Number of Studies",
       y = "Number of TF-target Pairs") +
  scale_x_continuous(breaks = 1:14) +  # Ensure x-axis has integer labels from 1 to 14
  theme_classic()

## 1.2 Summarize per TF ********************************************************************
### dot plot 
tf_summary <- tf_target_counts %>% 
  #filter(study_count >=5) %>%
  group_by(geneA) %>%
  summarize(
    num_reproducible_targets = sum(study_count > 5),  # Count targets appearing in more than 1 study
    total_targets = n(),
    proportion_reproducible = num_reproducible_targets / total_targets
  ) %>%
  arrange(desc(proportion_reproducible))



ggplot(tf_summary, aes(x = proportion_reproducible, y = reorder(geneA, proportion_reproducible))) +
#ggplot(tf_summary, aes(x = reorder(geneA, proportion_reproducible), y = proportion_reproducible)) +
  geom_point(color = "black", size = 2) +  # Dot plot with black points
  labs(
    title = "Proportion of Reproducible Targets per TF",
    x = "Proportion of geneB targets reproducible across studies",
    y = "TF (geneA)"
  ) +
  theme_classic() + 
  theme(
    axis.text.y = element_blank())
#So, if proportion_reproducible = 0.6, then 60% of the targets for that TF are observed 
# in at least two studies, making them reproducible 
# 20% of MEF2A’s geneB targets appear in 6 or more studies. 


# Create a frequency table and order by frequency
tf_freq <- as.data.frame(table(tf_target_counts$geneA)) %>%
  arrange(desc(Freq))  # Order TFs by descending frequency

# Plot
ggplot(tf_freq, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_point(color = "black", size = 2) +  # Dot plot
  labs(
    title = "Number of Targets per TF (geneA)",
    x = "TF (geneA)",
    y = "Frequency (Number of geneB targets)"
  ) + 
  theme_classic() + 
  theme(axis.text.x = element_blank())



#### update the dot plot above to color certain genes 
library(ggrepel)  # Load ggrepel for non-overlapping text labels

# Define the TFs to highlight
highlight_tfs <- c("FOXN3", "MEF2A", "RUNX1", "FOXP2", "JAZF1")

ggplot(tf_summary, aes(x = proportion_reproducible, y = reorder(geneA, proportion_reproducible))) +
  geom_point(aes(color = geneA %in% highlight_tfs), size = 2) +  # Color specific TFs
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +  # Red for selected TFs, black for others
  geom_text_repel(
    data = tf_summary %>% filter(geneA %in% highlight_tfs),  # Only label the highlighted TFs
    aes(label = geneA),
    color = "red",
    size = 4,  # Adjust label size
    nudge_x = 0.02,  # Small horizontal shift to avoid overlap
    direction = "both"
  ) +
  labs(
    title = "Proportion of Reproducible Targets per TF",
    x = "Proportion of geneB targets reproducible across studies",
    y = "TF (geneA)"
  ) +
  theme_minimal(base_size = 14) +  # Minimal theme
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white", color = NA),  # White background
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"  # Remove legend for simplicity
  )





## 1.3 cor dist ********************************************************************
ggplot(combined_df_top200, aes(x = rank_norm_corr, color = study)) +
  geom_density(alpha = 0.5, size = 0.5) +
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_classic() + 
  scale_color_viridis_d(option = "cividis") 

### facets 
ggplot(combined_df_top200, aes(x = rank_norm_corr)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~ study, scales = "free") +
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_minimal()


### TF-200 vs TF-allTrgs 
combined_df_top200$category <- "Top200"
combined_df_alltrg$category <- "AllTargets"
combined_df <- bind_rows(combined_df_top200, combined_df_alltrg)

test = combined_df %>% filter(study == "Mic_Hoffman-2023_CTL")

# Generate the density plot - one at a time (?)
ggplot(combined_df, aes(x = rank_norm_corr, color = category)) +
  geom_density(size = 0.5) +  # Set line size for clarity
  scale_color_manual(values = c("black", "firebrick")) + 
  facet_wrap(~ study, scales = "free") +  # Facet by study
labs(
  title = "Distribution of Rank Normalized Correlation Per Study",
  x = "Rank Normalized Correlation",
  y = "Density"
) +
  theme_bw() 


ggplot(test, aes(x = rank_norm_corr, color = category)) +
  geom_density(size = 0.5) +  # Set line size for clarity
  scale_color_manual(values = c("black", "firebrick"))
  #facet_wrap(~ study, scales = "free") +  # Facet by study
  labs(
    title = "Distribution of Rank Normalized Correlation Per Study",
    x = "Rank Normalized Correlation",
    y = "Density"
  ) +
  theme_bw() 


#  2. my analysis vs Alex's Mic =======================================================================================================================
alex_mic_allRank   = readRDS("/home/nelazzabi/rewiring/data/TFs/alex_mic_allRank.rds")
alex_mic_fishZ     = readRDS("/home/nelazzabi/rewiring/data/TFs/alex_mic_fishZ.rds")
  
sum(unique(combined_df_top200$geneA) %in% unique(alex_mic_allRank$geneA))
# [1] 848
sum(unique(combined_df_top200$geneA) %in% unique(alex_mic_fishZ$geneA))
# [1] 861


library(dplyr)
library(ggplot2)

# Sample datasets
# combined_df_top200 <- read.csv("path/to/combined_df_top200.csv")
# alex_mic_allRank <- read.csv("path/to/alex_mic_allRank.csv")

# Merge datasets on geneA and geneB
merged_df <- combined_df_top200 %>%
  inner_join(alex_mic_allRank, by = c("geneA", "geneB")) %>%
  rename(query_rank = rank_norm_corr, reference_rank = rank_correlation)

merged_dfFishZ <- combined_df_top200 %>%
  inner_join(alex_mic_fishZ, by = c("geneA", "geneB")) %>%
  rename(query_rank = rank_norm_corr, reference_rank = rank_correlation)

# Plot: Scatter plot per study
ggplot(merged_dfFishZ, aes(x = reference_rank, y = query_rank)) +
  geom_point(alpha = 0.1, color = "blue", size = 1.5) +
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "Comparison of Query vs. Reference (FZ) Rank Correlations",
    x = "Reference Rank Correlation",
    y = "Query Rank Normalized Correlation"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))

ggplot(merged_df, aes(x = reference_rank, y = query_rank)) +
  geom_point(alpha = 0.1, color = "blue", size = 1.5) +
  facet_wrap(~study, scales = "free") +  # Stratify by study
  theme_minimal() +
  labs(
    title = "Comparison of Query vs. Reference (allRank) Rank Correlations",
    x = "Reference Rank Correlation",
    y = "Query Rank Normalized Correlation"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))


ggplot(merged_df %>% filter(study == "Mic_Hoffman-2023_CTL"), aes(x = reference_rank, y = query_rank)) +
  geom_point(alpha = 0.6, color = "blue") +
  #geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  theme_minimal() +
  labs(
    title = "Comparison of Query vs. Reference (allRank) Rank Correlations",
    x = "Reference Rank Correlation",
    y = "Query Rank Normalized Correlation"
  ) +
  theme(strip.text = element_text(size = 12, face = "bold"))

 

#  3. my analysis & eRegulon of the 4 TFs =======================================================================================================================

my_data_TF_top200_AD_3TF <- list()
my_data_TF_top200_CTL_3TF <- list()

missing_in_control_per_study <- list()
missing_in_disease_per_study <- list()

controls_by_study <- list() 
disease_by_study <- list() 

for (i in 1:length(my_data_TF_allTrg_AD)) {
  study_AD  = names(my_data_TF_top200_AD)[i]
  study_CTL = names(my_data_TF_top200_CTL)[i]
  
  my_data_TF_top200_AD_3TF[[study_AD]]   <- my_data_TF_top200_AD[[study_AD]] %>% filter(geneA %in% gene_list$transcription_factors)   %>% mutate(pair = paste0(geneA, "_", geneB), condition = "Disease")
  my_data_TF_top200_CTL_3TF[[study_CTL]] <- my_data_TF_top200_CTL[[study_CTL]] %>% filter(geneA %in% gene_list$transcription_factors) %>% mutate(pair = paste0(geneA, "_", geneB), condition = "Control")
  
  missing_in_disease_per_study[[study_AD]]  <- setdiff(my_data_TF_top200_CTL_3TF[[study_CTL]]$pair, my_data_TF_top200_AD_3TF[[study_AD]]$pair)
  missing_in_control_per_study[[study_CTL]] <- setdiff(my_data_TF_top200_AD_3TF[[study_AD]]$pair, my_data_TF_top200_CTL_3TF[[study_CTL]]$pair)
  
  #print(length(missing_in_control_per_study[[study_CTL]]))
  #print(length(missing_in_disease_per_study[[study_AD]]))
  
  
  missing_in_disease_per_study[[study_AD]] <- as.data.frame(do.call(rbind, strsplit(missing_in_disease_per_study[[study_AD]], "_")))
  colnames(missing_in_disease_per_study[[study_AD]]) <- c("geneA", "geneB")
  
  missing_in_control_per_study[[study_CTL]] <- as.data.frame(do.call(rbind, strsplit(missing_in_control_per_study[[study_CTL]], "_")))
  colnames(missing_in_control_per_study[[study_CTL]]) <- c("geneA", "geneB")
  
  #print(dim(missing_in_control_per_study[[study_CTL]]))
  #print(dim(missing_in_disease_per_study[[study_AD]]))
  
  missing_control_data_study <- my_data_TF_allTrg_CTL[[study_CTL]] %>% filter(geneA %in% missing_in_control_per_study[[study_CTL]]$geneA & geneB %in% missing_in_control_per_study[[study_CTL]]$geneB)
  missing_disease_data_study <- my_data_TF_allTrg_AD[[study_AD]] %>%   filter(geneA %in% missing_in_disease_per_study[[study_AD]]$geneA & geneB %in% missing_in_disease_per_study[[study_AD]]$geneB)
  
  # Add the missing pairs to the control and disease datasets for the current study
  controls_by_study[[study_CTL]] <- bind_rows( my_data_TF_top200_CTL_3TF[[study_CTL]], missing_control_data_study)
  disease_by_study[[study_AD]]   <- bind_rows( my_data_TF_top200_AD_3TF[[study_AD]]  , missing_disease_data_study)
  
}



# Combine all studies back into one dataset
filtered_controls_alltrg <-  bind_rows(controls_by_study, .id = "study")
filtered_disease_alltrg <-   bind_rows(disease_by_study, .id = "study")

# Sort the datasets by study, geneA, and geneB to ensure order is consistent
filtered_controls_alltrg <- filtered_controls_alltrg[order(filtered_controls_alltrg$study, filtered_controls_alltrg$geneA, filtered_controls_alltrg$geneB), ]
filtered_disease_alltrg <- filtered_disease_alltrg[order(filtered_disease_alltrg$study, filtered_disease_alltrg$geneA, filtered_disease_alltrg$geneB), ]

# Check the final dimensions
dim(filtered_controls_alltrg)
dim(filtered_disease_alltrg)



filtered_combined <- bind_rows(filtered_controls_alltrg %>% mutate(condition = "Control"), filtered_disease_alltrg %>% mutate(condition = "Disease"))
filtered_combined$study = gsub("(_[^_]+)$", "", filtered_combined$study)

ggplot(filtered_combined, aes(x = rank_norm_corr, color = condition)) +
  geom_density(alpha = 0.5) + 
  scale_color_manual(values = c("black", "firebrick")) + 
  facet_wrap(~study, scales = "free") +  # Separate plots per study
  theme_minimal() +
  labs(title = "Distribution of Correlation Values per Study (Matching TF-Target Pairs)",
       x = "Rank Normalized Correlation",
       y = "Density",
       fill = "Condition") 





### ************************ original ************************
combined_df_alltrg_CTL = combined_df_alltrg
combined_df_alltrg_AD  = bind_rows(my_data_TF_allTrg_AD, .id = "study")

combined_df_top200_AD  <- bind_rows(my_data_TF_top200_AD, .id = "study")
combined_df_top200_CTL <- combined_df_top200

### 3TFs and all their targets that are in the eRegulon
filtered_controls_alltrg <- combined_df_top200_CTL %>% inner_join(eRegulonV1, by = c("geneA" = "TF", "geneB" = "Gene")) %>% mutate(condition = "Control")
filtered_disease_alltrg  <- combined_df_top200_AD  %>% inner_join(eRegulonV1, by = c("geneA" = "TF", "geneB" = "Gene")) %>% mutate(condition = "Disease")


filtered_controls_alltrgV1 = filtered_controls_alltrg %>% mutate(pair = paste0(geneA, "_", geneB))
filtered_disease_alltrgV1  = filtered_disease_alltrg  %>% mutate(pair = paste0(geneA, "_", geneB))

controls_by_study <- split(filtered_controls_alltrg, filtered_controls_alltrg$study)
disease_by_study  <- split(filtered_disease_alltrg, filtered_disease_alltrg$study)




# Split the datasets by study
controls_by_study <- split(filtered_controls_alltrg, filtered_controls_alltrg$study)
disease_by_study  <- split(filtered_disease_alltrg, filtered_disease_alltrg$study)

# Initialize lists to store missing data for each study
missing_in_control_per_study <- list()
missing_in_disease_per_study <- list()

# Loop through each study and find missing pairs
for(i in 1:length(controls_by_study)) {
  
  # Get control and disease pairs for the current study
  control_pairs_study <- unique(controls_by_study[[i]][, c("geneA", "geneB")])
  disease_pairs_study <- unique(disease_by_study[[i]][, c("geneA", "geneB")])
  
  # Get study names
  study_AD <- names(disease_by_study)[i]  # Correctly reference the study name
  study_CTL <- names(controls_by_study)[i]  # Correctly reference the study name
  
  # Find missing pairs for the current study
  missing_in_control_per_study[[study_CTL]] <- setdiff(disease_pairs_study, control_pairs_study)
  missing_in_disease_per_study[[study_AD]] <- setdiff(control_pairs_study, disease_pairs_study)
}

# Print out missing pairs for each study to review
missing_in_control_per_study
missing_in_disease_per_study

# Loop through the studies and add missing pairs to the datasets
for(i in 1:length(names(controls_by_study))) {  
  # Get the current study names
  study_AD <- names(disease_by_study)[i]
  study_CTL <- names(controls_by_study)[i]
  
  # Get the missing control and disease data for the current study
  missing_control_data_study <- combined_df_alltrg_CTL[combined_df_alltrg_CTL$study == study_CTL &
                                                         combined_df_alltrg_CTL$geneA %in% missing_in_control_per_study[[study_AD]]$geneA &
                                                         combined_df_alltrg_CTL$geneB %in% missing_in_control_per_study[[study_AD]]$geneB, ]
  
  missing_disease_data_study <- combined_df_alltrg_AD[combined_df_alltrg_AD$study == study_AD &
                                                        combined_df_alltrg_AD$geneA %in% missing_in_disease_per_study[[study_CTL]]$geneA &
                                                        combined_df_alltrg_AD$geneB %in% missing_in_disease_per_study[[study_CTL]]$geneB, ]
  
  # Add the missing pairs to the control and disease datasets for the current study
  controls_by_study[[study_CTL]] <- bind_rows(controls_by_study[[study_CTL]], missing_control_data_study)
  disease_by_study[[study_AD]] <- bind_rows(disease_by_study[[study_AD]], missing_disease_data_study)
}

# Combine all studies back into one dataset
filtered_controls_alltrg <- bind_rows(controls_by_study)
filtered_disease_alltrg <- bind_rows(disease_by_study)

# Sort the datasets by study, geneA, and geneB to ensure order is consistent
filtered_controls_alltrg <- filtered_controls_alltrg[order(filtered_controls_alltrg$study, filtered_controls_alltrg$geneA, filtered_controls_alltrg$geneB), ]
filtered_disease_alltrg <- filtered_disease_alltrg[order(filtered_disease_alltrg$study, filtered_disease_alltrg$geneA, filtered_disease_alltrg$geneB), ]

# Check the final dimensions
dim(filtered_controls_alltrg)
dim(filtered_disease_alltrg)



filtered_combined <- bind_rows(filtered_controls_alltrg %>% mutate(condition = "Control"), filtered_disease_alltrg %>% mutate(condition = "Disease"))
filtered_combined$study = gsub("(_[^_]+)$", "", filtered_combined$study)

ggplot(filtered_combined, aes(x = rank_norm_corr, color = condition)) +
  geom_density(alpha = 0.5) + 
  scale_color_manual(values = c("black", "firebrick")) + 
  facet_wrap(~study, scales = "free") +  # Separate plots per study
  theme_minimal() +
  labs(title = "Distribution of Correlation Values per Study (Matching TF-Target Pairs)",
       x = "Rank Normalized Correlation",
       y = "Density",
       fill = "Condition") 




























filtered_controls_alltrg <- combined_df_alltrg_CTL %>% inner_join(eRegulonV1, by = c("geneA" = "TF", "geneB" = "Gene")) %>% mutate(condition = "Control")
filtered_disease_alltrg  <- combined_df_alltrg_AD  %>% inner_join(eRegulonV1, by = c("geneA" = "TF", "geneB" = "Gene")) %>% mutate(condition = "Disease")

all(unique(filtered_controls_alltrg$geneA) %in% unique(filtered_disease_alltrg$geneA))
#[1] TRUE
all(unique(filtered_controls_alltrg$geneB) %in% unique(filtered_disease_alltrg$geneB))
#[1] TRUE

filtered_combined <- bind_rows(filtered_controls_alltrg %>% mutate(condition = "Control"), filtered_disease_alltrg %>% mutate(condition = "Disease"))
filtered_combined$study = gsub("(_[^_]+)$", "", filtered_combined$study)

ggplot(filtered_combined, aes(x = rank_norm_corr, color = condition)) +
  geom_density(alpha = 0.5) + 
  scale_color_manual(values = c("black", "firebrick")) + 
  facet_wrap(~study, scales = "free") +  # Separate plots per study
  theme_minimal() +
  labs(title = "Distribution of Correlation Values per Study (Matching TF-Target Pairs)",
       x = "Rank Normalized Correlation",
       y = "Density",
       fill = "Condition") 




### work with the top and then create a list of all possible gene-pairs and then extract their corr from the allTrg files 

  
  
  
  
  
  
  
  
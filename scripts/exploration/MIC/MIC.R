

##### read in data ========================================================================================================================
##### =====================================================================================================================================
library(data.table)

### interesting genes 
TFs = readRDS("/home/nelazzabi/rewiring/data/TFs/02_Human_TR_coexpr_rankings.RDS") 
ADgene_list <- c(
  "ABCA7", "ABI3", "ACE2", "ADAM10", "APH1B", "APOE", "BIN1", "CD2AP", 
  "CR1", "ECHDC3", "FERMT2", "HS3ST1", "IL34", "MINK1", "PLCG2", "PLD3", 
  "PRKD3", "PTK2B", "SCIMP", "SHARPIN", "SLC2A4", "SORL1", "TREM2", "UNC5C"
)

input_files <- list.files(path = "/cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list", pattern = "Mic_.*_edgeList.csv", full.names = TRUE)
input_files <- list.files(path = "/cosmos/data/project-data/NW-rewiring/coExpr/disease/edge-list", pattern = "Mic_.*_corEdgeLists.csv", full.names = TRUE)


test_TF_allTargets = process_file(input_files[2], TFs, top_n_targets = NULL) 
test_TF_topTargets = process_file(input_files[2], TFs, top_n_targets = 200) 



View(table(test_TF_allTargets$geneA))
View(table(test_TF_topTargets$geneA))





df = fread(input_files[1])

df2 = df[, pair := paste0(pmin(geneA, geneB), "_", pmax(geneA, geneB))]

df_unique <- unique(df2, by = "pair")
df_unique[, pair := NULL]



View(table(testTest1$geneA))
















tf_top200_list <- setNames(
  lapply(input_files, function(file) extract_TF_top200(file, names(TFs))),
  gsub("_edgeList\\.csv$", "", basename(input_files))
)








d_list_top200 <- setNames(lapply(input_files, process_file), gsub("_edgeList\\.csv$", "", basename(input_files)))
d_list_alltrg <- setNames(lapply(input_files, process_fileAll), gsub("_edgeList\\.csv$", "", basename(input_files)))


# sanity check 1 = there are same TF in all studies 
# # of targets per TF is the same 
# dim of each study is the same 

top200_AD  <- d_list_top200[grepl("_AD", names(d_list_top200))]
top200_CTL <- d_list_top200[grepl("_CTL", names(d_list_top200))]

alltrg_AD  <- d_list_alltrg[grepl("_AD", names(d_list_alltrg))]
alltrg_CTL <- d_list_alltrg[grepl("_CTL", names(d_list_alltrg))]








### sanity check 
### how many TFs are present in my analysis 
for (i in 1:length(filtered_data_CTL)) {
  #print(length(unique(filtered_data_CTL[[i]]$geneA)))
  #print(length(unique(filtered_data_CTL[[i]]$geneB))) #this is not really that important 
  print(dim(filtered_data_CTL[[i]]))
}

View(table(TF_allTarget_CTL[[2]]$geneB))








##### CTL summarize =======================================================================================================================
##### =====================================================================================================================================

library(dplyr)

### wranlge data 
# 1. Combine data into one dataframe
combined_df <- bind_rows(filtered_data_CTL, .id = "study")
combined_df_alltr <- bind_rows(TF_allTarget_CTL, .id = "study")

combined_df100 <- bind_rows(filtered_data_CTL_100, .id = "study")



# 2. Count occurrences of each TF-target pair
tf_target_counts <- combined_df %>%
  group_by(geneA, geneB) %>%
  summarise(freq = n(), .groups = "drop")

tf_target_counts100 <- combined_df100 %>%
  group_by(geneA, geneB) %>%
  summarise(freq = n(), .groups = "drop")

# 3. Count how many TF-target pairs appear in exactly X studies - per TF 
tf_summary <- tf_target_counts %>%
  group_by(geneA) %>%
  count(freq) %>%
  pivot_wider(names_from = freq, values_from = n, values_fill = 0)

tf_summary100 <- tf_target_counts100 %>%
  group_by(geneA) %>%
  count(freq) %>%
  pivot_wider(names_from = freq, values_from = n, values_fill = 0)




filtered_tf_target_counts <- tf_target_counts %>%
  filter(freq >= 5)  # Keep only pairs appearing in 5 or more studies

filtered_tf_target_counts100 <- tf_target_counts100 %>%
  filter(freq >= 5)  # Keep only pairs appearing in 5 or more studies

tf_summaryV1 <- filtered_tf_target_counts %>%
  group_by(geneA) %>%
  count(freq) %>%
  pivot_wider(names_from = freq, values_from = n, values_fill = 0)

tf_summaryV1_100 <- filtered_tf_target_counts100 %>%
  group_by(geneA) %>%
  count(freq) %>%
  pivot_wider(names_from = freq, values_from = n, values_fill = 0)


# Convert to matrix format
filtered_tf_matrix <- filtered_tf_target_counts %>%
  tidyr::pivot_wider(names_from = geneB, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames("geneA") %>%
  as.matrix()

filtered_tf_matrix100 <- filtered_tf_target_counts100 %>%
  tidyr::pivot_wider(names_from = geneB, values_from = freq, values_fill = 0) %>%
  tibble::column_to_rownames("geneA") %>%
  as.matrix()



##### CTL visualize =======================================================================================================================
##### =====================================================================================================================================

library(ggplot2)

### 1. bar plots 
### top 200 co-expressed partners 
ggplot(tf_target_counts, aes(x = freq)) +
  geom_bar(fill = "blue", alpha = 0.8) +
  labs(title = "Distribution of TF-Target Pairs Across Studies",
       x = "Number of Studies",
       y = "Number of TF-Target Pairs") +
  scale_x_continuous(breaks = 1:14) +  # Ensure x-axis has integer labels from 1 to 14
  theme_minimal()

ggplot(filtered_tf_target_counts, aes(x = freq)) +
  geom_bar(fill = "blue", alpha = 0.8) +
  labs(title = "Distribution of TF-Target Pairs Across Studies",
       x = "Number of Studies",
       y = "Number of TF-Target Pairs") +
  scale_x_continuous(breaks = 1:14) +  # Ensure x-axis has integer labels from 1 to 14
  theme_minimal()


### top 100 co-expressed partners 
ggplot(tf_target_counts100, aes(x = freq)) +
  geom_bar(fill = "blue", alpha = 0.8) +
  labs(title = "Distribution of TF-Target Pairs Across Studies",
       x = "Number of Studies",
       y = "Number of TF-Target Pairs") +
  scale_x_continuous(breaks = 1:14) +  # Ensure x-axis has integer labels from 1 to 14
  theme_minimal()

ggplot(filtered_tf_target_counts100, aes(x = freq)) +
  geom_bar(fill = "blue", alpha = 0.8) +
  labs(title = "Distribution of TF-Target Pairs Across Studies",
       x = "Number of Studies",
       y = "Number of TF-Target Pairs") +
  scale_x_continuous(breaks = 1:14) +  # Ensure x-axis has integer labels from 1 to 14
  theme_minimal()


### 2. heatmap of TFs vs study count frequency 
pheatmap::pheatmap(as.matrix(tf_summary[,-1]), 
                   cluster_rows = TRUE, 
                   cluster_cols = FALSE, 
                   color = viridis::viridis(100, option = "inferno"), 
                   #color = viridis::viridis(100, option = "viridis", direction = -1), 
                   main = "top 200 TF-Target Pair Frequency Across Studies")


pheatmap::pheatmap(as.matrix(tf_summaryV1[,-1]), 
                   cluster_rows = TRUE, 
                   cluster_cols = FALSE, 
                   color = viridis::viridis(100, option = "inferno"), 
                   #color = viridis::viridis(100, option = "viridis", direction = -1), 
                   main = "top 200 TF-Target Pair Frequency Across Studies")



### top 100 co-expressed partners 
pheatmap::pheatmap(as.matrix(tf_summary100[,-1]), 
                   cluster_rows = TRUE, 
                   cluster_cols = FALSE, 
                   color = viridis::viridis(100, option = "inferno"), 
                   #color = viridis::viridis(100, option = "viridis", direction = -1), 
                   main = "top 100 TF-Target Pair Frequency Across Studies")


pheatmap::pheatmap(as.matrix(tf_summaryV1_100[,-1]), 
                   cluster_rows = TRUE, 
                   cluster_cols = FALSE, 
                   color = viridis::viridis(100, option = "inferno"), 
                   #color = viridis::viridis(100, option = "viridis", direction = -1), 
                   main = "top 100 TF-Target Pair Frequency Across Studies")


### 3. heatmap of 5> occurant TF-targets , Tfs vs Targets 
library(viridis)
pheatmap::pheatmap(filtered_tf_matrix,
         cluster_rows = TRUE, cluster_cols = TRUE,
         color = viridis::viridis(100),
         main = "TF-Target Pair Frequency (≥30% of Studies)")

pheatmap::pheatmap(filtered_tf_matrix100,
                   cluster_rows = TRUE, cluster_cols = TRUE,
                   color = viridis::viridis(500, option = "inferno"),
                   main = "TF-Target Pair Frequency (≥30% of Studies)")


### 4. density plot per TF 
library(ggplot2)
ggplot(filtered_tf_target_counts, aes(x = freq, color = geneA)) +
  geom_density(alpha = 0.6, size = 0.6) +  
  labs(title = "Density of TF-Target Pair Frequencies per TF",
       x = "Number of Studies",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:14) + 
  theme(legend.position = "none")  # Hide legend

ggplot(filtered_tf_target_counts100, aes(x = freq, color = geneA)) +
  geom_density(alpha = 0.6, size = 0.6) +  
  labs(title = "Density of TF-Target Pair Frequencies per TF",
       x = "Number of Studies",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:14) + 
  theme(legend.position = "none")  # Hide legend



##### extract frequent pairs ==============================================================================================================
##### =====================================================================================================================================

#### extract the most reproducible co-expressed targets per TF 
# 1. filter out TF-target pairs that appear in at least 5 studies 
filtered_tf_target_counts 

# what is the total #of targets per TF 
# what is the table summary of targets (are they shared among TFs?) 

# 1. Table: Total Number of Targets per TF (geneA)
tf_summary <- filtered_tf_target_counts %>%
  group_by(geneA) %>%
  summarise(total_targets = n(), .groups = "drop")

# 2. Table: Frequency of Each Target (geneB) Across TFs
target_summary <- filtered_tf_target_counts %>%
  group_by(geneB) %>%
  summarise(occurrence = n(), .groups = "drop")

hist(tf_summary$total_targets, breaks = 40)
hist(target_summary$occurrence, breaks = 60)





##### cor distribution for the frequent pairs ============================================================================================
##### =====================================================================================================================================
library(dplyr)
library(ggplot2)
library(purrr)


extracted_correlation <- map(filtered_data_CTL, function(study_df) {
  study_df %>%
    inner_join(filtered_tf_target_counts, by = c("geneA", "geneB")) %>%
    select(geneA, geneB, rank_normalized_correlation)
})
combined_correlation <- bind_rows(extracted_correlation, .id = "Study")
combined_correlation_filt <- combined_correlation %>% filter(rank_normalized_correlation > 0.8)
#tf_target_counts_freq <- combined_correlation_filt %>%
  #group_by(geneA, geneB) %>%
  #summarise(freq = n(), .groups = "drop")


length(unique(combined_correlation_filt$geneA))
length(unique(combined_correlation_filt$geneB))

### all in one plot 
ggplot(combined_correlation_filt, aes(x = rank_normalized_correlation, fill = Study)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_minimal() #+
  #theme(legend.position = "none")  # Remove legend if needed

# OR Extended color-blind friendly palette with 14 colors ******
cb_friendly_palette_extended <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                                  "#D55E00", "#CC79A7", "#F1C6D1", "#B3B3B3", "#EAC300",
                                  "#A0D3E5", "#8E00C6", "#91BF5E", "#8E44AD")

ggplot(combined_correlation, aes(x = rank_normalized_correlation, fill = Study)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = cb_friendly_palette_extended) +  # Assign distinct colors to each study
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_minimal()


### facets 
ggplot(combined_correlation_filt, aes(x = rank_normalized_correlation)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~ Study, scales = "free") +
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_minimal()









##### compare frequent pairs to Alex's global signal ======================================================================================
##### =====================================================================================================================================
library(dplyr)
library(purrr)

filtered_tf_target_counts # all cor 

combined_correlation_filt # only 0.8 cor
 
filtered_tf_target_countsTop <- filtered_tf_target_counts %>%
  semi_join(combined_correlation_filt, by = c("geneA", "geneB"))

# Add Top200_count column
filtered_tf_target_counts_integrated <- filtered_tf_target_countsTop %>%
  rowwise() %>%
  mutate(Top200_count = if_else(
    geneA %in% names(TFs) & geneB %in% rownames(TFs[[geneA]]),
    TFs[[geneA]][geneB, "Top200_count"],  # Extract Top200_count
    NA_real_  # Assign NA if not found
  )) %>%
  ungroup()


count_zeros <- sum(filtered_tf_target_counts_integrated$Top200_count == 0, na.rm = TRUE)
count_NAs <- sum(is.na(filtered_tf_target_counts_integrated$Top200_count))

long_data <- filtered_tf_target_counts_integrated %>%
  pivot_longer(cols = c(freq, Top200_count),
               names_to = "Measure",
               values_to = "Value")

### visualize ******
library(dplyr)
library(tidyr)
library(ggplot2)

ggplot(long_data %>% filter(!is.na(Value) & Value != 0), aes(x = Value, fill = Measure)) +
  geom_density(alpha = 0.5) +  # alpha for transparency
  scale_fill_manual(values = c("skyblue", "orange")) +  # Custom colors
  labs(title = "Density Plot of freq and Top200_count",
       x = "Value",
       y = "Density",
       fill = "Measure") +
  theme_minimal()

# 1. what is the Separate freq distribution for each measure

filtered_data <- long_data %>%
  filter(!is.na(Value) & Value != 0)

# Plot for 'freq'
ggplot(filtered_data %>% filter(Measure == "freq"), aes(x = Value, fill = Measure)) +
  geom_histogram(binwidth = 1, alpha = 0.5, position = "identity", color = "black", size = 0.1) +
  scale_fill_manual(values = c("skyblue")) +  # Custom color for 'freq'
  labs(title = "Histogram of freq",
       x = "Value",
       y = "Count",
       fill = "Measure") +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold"))

# Plot for 'Top200_count'
ggplot(filtered_data %>% filter(Measure == "Top200_count"), aes(x = Value, fill = Measure)) +
  geom_histogram(binwidth = 1, alpha = 0.5, position = "identity", color = "black", size = 0.1) +
  scale_fill_manual(values = c("orange")) +  # Custom color for 'Top200_count'
  labs(title = "Histogram of Top200_count",
       x = "Value",
       y = "Count",
       fill = "Measure") +
  theme_minimal() +
  theme(axis.text = element_text(face = "bold"))


# 2. what is the freq distriution for the ones that have 0s in Alex's work 
# 3. what is the freq distribution for the ones that have NAs in Alex's work 
zero_pairs <- filtered_tf_target_counts_integrated %>% filter(Top200_count == 0)
na_pairs <- filtered_tf_target_counts_integrated %>% filter(is.na(Top200_count))

# Extract correlations from 'combined_correlation_filt'
zero_pairs_corr <- combined_correlation_filt %>% filter(paste(geneA, geneB) %in% paste(zero_pairs$geneA, zero_pairs$geneB))
na_pairs_corr <- combined_correlation_filt %>% filter(paste(geneA, geneB) %in% paste(na_pairs$geneA, na_pairs$geneB))

### OR 
zero_pairs_corr <- combined_correlation %>% filter(paste(geneA, geneB) %in% paste(zero_pairs$geneA, zero_pairs$geneB))
na_pairs_corr <- combined_correlation %>% filter(paste(geneA, geneB) %in% paste(na_pairs$geneA, na_pairs$geneB))

ggplot(zero_pairs_corr, aes(x = rank_normalized_correlation)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~ Study, scales = "free") +
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_minimal()

# 4. what is the dist of the number of unique pairs per study 
library(ggplot2)
library(dplyr)
library(tidyr)

# Combine the two dataframes
combined_pairs_corr <- bind_rows(
  zero_pairs_corr %>% mutate(Type = "Zero Top200 Count"),
  na_pairs_corr %>% mutate(Type = "NA Top200 Count")
)

# Count the number of unique TF-target pairs per study
unique_pairs_per_study <- combined_pairs_corr %>%
  group_by(Study, Type) %>%
  summarise(unique_pairs = n_distinct(paste(geneA, geneB)), .groups = "drop")

# Reshape data to make it easier to visualize
unique_pairs_wide <- unique_pairs_per_study %>%
  pivot_wider(names_from = Type, values_from = unique_pairs, values_fill = list(unique_pairs = 0))

# Sort the data by 'Zero Top200 Count' in descending order
unique_pairs_wide <- unique_pairs_wide %>%
  arrange(desc(`Zero Top200 Count`)) 

# Convert 'Study' to a factor, ordered by 'Zero Top200 Count'
unique_pairs_wide$Study <- factor(unique_pairs_wide$Study, levels = unique_pairs_wide$Study)

# Plot the result
ggplot(unique_pairs_wide, aes(x = Study)) +
  geom_bar(aes(y = `Zero Top200 Count`, fill = "Zero Top200 Count"), stat = "identity", position = "dodge") +
  geom_bar(aes(y = `NA Top200 Count`, fill = "NA Top200 Count"), stat = "identity", position = "dodge") +
  labs(title = "Number of Unique TF-Target Pairs per Study",
       x = "Study",
       y = "Unique TF-Target Pairs") +
  scale_fill_manual(values = c("Zero Top200 Count" = "skyblue", "NA Top200 Count" = "orange")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(face = "bold")) + 
  scale_y_continuous(expand = c(0, 0))  # Remove extra space on y-axis +














##### checking AD genes  ==================================================================================================================
##### =====================================================================================================================================
library(dplyr)


# Count AD genes in geneA for each dataset
ad_gene_counts <- purrr::map_dfr(names(filtered_data_CTL), function(study) {
  data <- filtered_data_CTL[[study]]
  
  # Identify AD genes in geneA and geneB
  ad_geneA_overlap <- unique(data$geneA[data$geneA %in% ADgene_list])
  ad_geneB_overlap <- unique(data$geneB[data$geneB %in% ADgene_list])
  
  # Count occurrences
  count_geneA <- length(ad_geneA_overlap)
  count_geneB <- length(ad_geneB_overlap)
  count_geneB_Perc <- count_geneB / length(ADgene_list)
  
  # Create tibble with new column listing AD genes in geneB
  tibble(
    Study = study, 
    AD_GeneA_Count = count_geneA, 
    AD_GeneB_Count = count_geneB, 
    AD_GeneB_Count_Perc = count_geneB_Perc,
    AD_Genes_in_GeneB = paste(ad_geneB_overlap, collapse = ", ")
  )
})

ad_gene_counts %>% arrange(desc(AD_GeneB_Count))



# count times AD gene in GeneB across studies
library(tidyr)
ad_geneB_summary <- ad_gene_counts %>%
  filter(AD_Genes_in_GeneB != "") %>%  # Remove empty entries
  separate_rows(AD_Genes_in_GeneB, sep = ", ") %>%  # Split multiple genes into separate rows
  count(AD_Genes_in_GeneB, name = "Study_Count") %>%  # Count occurrences across studies
  arrange(desc(Study_Count))  # Sort by frequency



### extract cor values for those targets across studies 
AD_geneB_cor <- combined_correlation %>%
  filter(geneB %in% ad_geneB_summary$AD_Genes_in_GeneB)

### plot it on top of the overall cor distribution (prior to doing any filtering) 
ggplot(AD_geneB_cor, aes(x = rank_normalized_correlation)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(~ Study, scales = "free") +
  labs(title = "Distribution of Rank Normalized Correlation Per Study",
       x = "Rank Normalized Correlation",
       y = "Density") +
  theme_minimal()






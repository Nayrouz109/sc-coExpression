


df <- fread("/cosmos/data/project-data/NW-rewiring/data/unibind/Curated_targets_all_Sept2023.tsv")

### harmonize data across species ==================================================================
library(dplyr)
library(stringr)

# 1. Harmonize gene symbols to uppercase (as a proxy for human-style symbols)
df_clean <- df %>%
  mutate(
    TF_Symbol_std = toupper(TF_Symbol),
    Target_Symbol_std = toupper(Target_Symbol)
  )

# 2. Collapse species-specific rows into unique TF–target pairs
df_collapsed <- df_clean %>%
  group_by(TF_Symbol_std, Target_Symbol_std) %>%
  summarise(
    TF_Species = paste(unique(TF_Species), collapse = ";"),
    Target_Species = paste(unique(Target_Species), collapse = ";"),
    Cell_Type = paste(unique(Cell_Type), collapse = ";"),
    Experiment_Type = paste(unique(Experiment_Type), collapse = ";"),
    PubMed_ID = paste(unique(PubMed_ID), collapse = ";"),
    Databases = paste(unique(Databases), collapse = ";"),
    .groups = "drop"
  )

# 3. Rename columns for clarity
df_ref <- df_collapsed %>%
  rename(
    TF_Symbol = TF_Symbol_std,
    Target_Symbol = Target_Symbol_std
  )



saveRDS(df_ref, "/home/nelazzabi/rewiring/data/modalities/lit_curated_targets.rds")



##################    testing   ########################################################################3

intersect(
  unique(df_ref %>% filter(TF_Symbol == "RUNX1") %>% pull(Target_Symbol)),
  tf_top_with_binding %>% filter(geneA == "RUNX1") %>% pull(geneB)
)


length(unique((df %>% filter(TF_Species == "Human"))$TF_Symbol))




test = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/Ast_t200_rpr_tf.rds")

data_path <- "/home/nelazzabi/rewiring/data/coExpr/0all-human"
files <- list.files(path = data_path, pattern = "t200_.*tf\\.rds$", full.names = TRUE)


read_and_format <- function(file) {
  df <- readRDS(file) %>%
    select(geneA, meanTop200Ctl, mic_TF, brain_TF)
  
  # Extract cell type from filename (e.g., "Ast" from "Ast_t200_rpr_tf.rds")
  cell_type <- str_extract(basename(file), "^[^_]+")
  mean_col <- paste0("meanTop200Ctl_", cell_type)
  
  # Rename column
  df <- df %>%
    rename(!!mean_col := meanTop200Ctl)
  
  return(df)
}


merged_df <- reduce(
  lapply(files, read_and_format),
  full_join,
  by = c("geneA", "mic_TF", "brain_TF")
)

# View result
print(merged_df)


saveRDS(merged_df, "/home/nelazzabi/rewiring/data/coExpr/0all-human/intgr/TR_rpr_summary.rds")










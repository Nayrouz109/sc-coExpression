

library(data.table)


test = readRDS("/space/scratch/singlecell_DEA_jan2025/sc_data_compiled_feb2025_v2.rds")


meta_dt <- test$mouse$meta

# Initialize an empty list to collect subsetFactor entries
subset_summary <- rbindlist(lapply(1:nrow(meta_dt), function(i) {
  exp_id <- meta_dt$experiment.ID.x[i]
  subset_dt <- meta_dt$subsetFactor[[i]]
  
  # Add the experiment.ID.x as a column to each subsetFactor entry
  subset_dt[, experiment.ID.x := exp_id]
  
  return(subset_dt)
}), fill = TRUE)

# Now, optionally, summarize the subset information per experiment.ID.x
summary_dt <- subset_summary[, .(
  subset_factor_terms = paste(unique(factor.category), collapse = "; "),
  subset_factor_values = paste(unique(value), collapse = "; ")
), by = experiment.ID.x]

# Long-format summary: one row per experiment.ID.x and subsetFactor entry
summary_dt <- subset_summary[, .(
  experiment.ID.x,
  subset_factor_term = factor.category,   # or use correct column name
  subset_factor_value = value      # or use correct column name
)] %>% unique() 


# View the result
head(summary_dt)

# Gene coExpression "Network" Rewiring

### Project purpose and goals 
This project curates large-scale human brain sc RNA-seq datasets and computes **cell-type-specific coexpression networks**. It then compares networks to:

- Study patterns of coexpression that are **cell type-specific**
- Study patterns that are **shared across cell types**
- Study gene-gene coexpression relationships that may **change in Alzheimer’s disease**

### Project summary
The project develops a framework for large-scale sc RNA-seq integration spanning **>20M cells across cohorts**. It:

- Engineers a scalable Python pipeline for correlation-based gene coexpression analysis, using rank-based and Spearman methods tailored to sparse scRNA-seq data
- Builds an end-to-end workflow for signal reproducibility benchmarking and cross-species validation
- Applies cross-validation and resampling methods to assess coexpression signal robustness
- Uses network-based analytics to map gene-module dysregulation in Alzheimer’s disease

---
### Where to find data and scripts 
- **data**:   /cosmos/data/project-data/NW-rewiring/data
- **coExpr matrices**: /cosmos/data/project-data/NW-rewiring/coExpr
- **script**: /home/nelazzabi/rewiring/scripts


## Project Structure

This repository contains the **`/scripts` folder**, which includes all scripts needed to run the analyses. The scripts are organized in different subfolders for clarity:

### `scripts/analysis`
Contains scripts that perform all analyses:

- **`scripts/analysis/1.only-cor`**  
  Computes correlation matrices using user-defined parameters and different analysis modes (across cells, across subjects, etc.)
  
- **`scripts/analysis/3.topCoExp/pipeline`**  
  Performs signal reproducibility benchmarking, cross-validation, and resampling analyses

- **`scripts/analysis/8.plusModalities`**  
  Integrates additional data modalities (e.g., binding evidence) into coexpression analyses

### `scripts/data-wrangling`
Contains scripts for wrangling large-scale human brain single-cell datasets into a unified format for downstream analysis.

### `scripts/functions`
Contains utility functions that are sourced in scripts across the `analysis` folder.

### Other Folders
- **`exploration`**: miscellaneous exploration scripts  
- **`util`**: general helper scripts  
- **`vis-summarize`**: visualization and summarization scripts  

---

### Contact 
for questions, contact nayrouz@student.ubc.ca


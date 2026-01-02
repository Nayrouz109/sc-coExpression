##### SUBJECT/SAMPLE LEVEL COEXPRESSION  #######################################
################################################################################
## Input: Pseudobulked samples in .csv format
################################################################################
##### USAGE  ###################################################################
#### RScript /home/nelazzabi/rewiring/scripts/analysis/7.cellvsample/xSample_coex.R <script>
#### options : xcell_reprod_prep 
############## xsubj_reprod_prep
############## compute_fisher_xcell
############## compute_fisher_xsubj


  library(tidyverse)
  
  ##### EDIT THIS PART TO SELECT DATASETS/FILES
  dataset <- c('ling-2024','ROSMAP-Mathys-2024','MSBB-2024', 'RADC-2024','HBCC-2024','SZBDMulti-2024',
               'MSSM-2024','Rexach-2024','ROSMAP-erosion-2023','Kamath-NPH-2023',
               'Hoffman-2023','ROSMAP-Mathys-2019','ROSMAP-Mathys-2023','ROSMAP-TREM2-2020',
               'ROSMAP-DeJager-2024','ROSMAP-DeJager-2023','SEA-AD-MTG-2024','SEA-AD-DLPFC-2024')
  cell_types <- c('Mic', 'Oli', 'Opc', 'Exc', 'Inh', 'Ast')
  
  
  ##### INPUT SHOULD BE DIR TO PSEUDOBULKED SAMPLES IN '.csv' FORMAT
  input_dir <- '/space/scratch/bxu_project/data/pseudobulk_cpm/control/'
  
  ##### EDIT OUTPUT DIR 
  output_dir <- '/space/scratch/bxu_project/pavlidislab_project/pan/subj_coex/redo'
  
  
  
  computeCoexMat <- function(exprmat) {
    exprmat <- exprmat[order(rownames(exprmat)), order(colnames(exprmat))]
    coexMat <- cor(t(exprmat), method = "pearson")
    return(coexMat)
  }
  
  normalizeCoexmat <- function(coexmat) {
    
    rankCoexmat <- matrix(nrow = nrow(coexmat), ncol = ncol(coexmat))
    rownames(rankCoexmat) <- rownames(coexmat)
    colnames(rankCoexmat) <- colnames(coexmat)
    
    lowerCoords <- which(lower.tri(rankCoexmat), arr.ind = TRUE) # use the lower tri indices as primary
    upperCoords <- lowerCoords[, c("col", "row")] # figure out the matching coordinates to "flip it back"
    colnames(upperCoords) <- c("row", "col")
    
    currCoexVec <- coexmat[lowerCoords]
    currCoexVec[is.na(currCoexVec)] <- currCoexVec %>% median(na.rm = TRUE)
    
    coexRanks <- currCoexVec %>% rank(ties.method = "average")
    coexRanks <- coexRanks / max(coexRanks)
    
    rankCoexmat[lowerCoords] <- coexRanks
    rankCoexmat[upperCoords] <- coexRanks
    diag(rankCoexmat) <- 1
    
    return(rankCoexmat)
  }
  
  computeCoexMats <- function(datasets) {
    
    allfiles <- list.files(input_dir, full.names = TRUE)
    # allfiles <- allfiles[grepl("expr", basename(allfiles))]
    
    # for every dataset in list
    for (d in datasets) {
      for (cell in cell_types) {
        exprfiles <- allfiles[grepl(d, basename(allfiles)) & grepl(cell, basename(allfiles))]
        if (length(exprfiles)==0) {
          next
        }
        for (f in exprfiles) {
          print(paste("Processing file:", f))
          expr <- read_csv(f)
          expr <- expr %>% column_to_rownames(var = colnames(expr)[1])
          
          
          expr <- expr[order(rownames(expr)), ]
          coex <- computeCoexMat(expr)
          coex <- normalizeCoexmat(coex)
          out_file <- file.path(output_dir, paste0(cell,'_', d, "_sbj_coexmats", ".csv"))
          write_csv(as.data.frame(coex), out_file)
          print(paste("Saved coexpression matrix to", out_file))
          
        }
      }
    }
  }
  
  computeCoexMats(dataset)





args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript file_name.R <function_name>\nAvailable functions: ")
}

script_name <- args[1]

# Run the requested function if it exists
if (exists(script_name) && is.function(get(script_name))) {
  do.call(script_name, list())  
} else {
  stop(paste("Error: Invalid function name:", script_name))
}

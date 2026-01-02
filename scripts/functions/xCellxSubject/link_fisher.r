

# ============================================================================================================================================================
# ============================================================================================================================================================
# ============================================================================================================================================================
# LINK TABLE  
# ============================================================================================================================================================
# ============================================================================================================================================================
# ============================================================================================================================================================


xCellxSubj_reprod_prep <- function() {
  cat("Running xCellxSubj_reprod_prep...")
  
  library(tidyverse)
  library(parallel)
  
  dataset <- c('Ling-2024','ROSMAP-Mathys-2024','MSBB-2024', 'RADC-2024','HBCC-2024','SZBDMulti-2024',
               'MSSM-2024','Rexach-2024','ROSMAP-erosion-2023','Kamath-NPH-2023',
               'Hoffman-2023','ROSMAP-Mathys-2019','ROSMAP-Mathys-2023','ROSMAP-TREM2-2020',
               'ROSMAP-DeJager-2024','ROSMAP-DeJager-2023','SEA-AD-MTG-2024','SEA-AD-DLPFC-2024')
  cell_types <- c('Mic', 'Oli', 'Opc', 'Exc', 'Inh', 'Ast')
  
  # input directory containing xcell coexpression matrices per dataset
  input_dir <- '/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/edge-list/all-matrix'
  output_dir <- '/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/lnkTables'
  
  
  coexmats <- list.files(input_dir, full.names = TRUE)
  
  
  getPairIds <- function(coexMat) {
    geneA <- rownames(coexMat)[row(coexMat)[lower.tri(coexMat)]]
    geneB <- colnames(coexMat)[col(coexMat)[lower.tri(coexMat)]]
    pairId <- paste0(geneA, ".", geneB)
    return(pairId)
  }
  
  vectorize <- function(coexMat) {
    coexMat[lower.tri(coexMat, diag = FALSE)]
  }
  
  
  for (d in dataset) {
    for (cell in cell_types) {
      cat('Processing:', d, cell, '\n')
      mats <- coexmats[grepl(d, basename(coexmats), ignore.case = TRUE) & grepl(cell, basename(coexmats), ignore.case = TRUE)]
      
      
      lapply(mats, function(lnkmat) {
        dataset_name <- d
        
        lnkmat <- tryCatch(read_rds(lnkmat), 
                           error = function(e) { cat("Error reading:", lnkmat, "\n"); return(NULL) })
        if (is.null(lnkmat)) return(NULL)  # Skip on read error
             
        lnkmat[lnkmat >= 0.99] <- 1
        lnkmat[lnkmat < 0.99] <- 0
        
        lnktbl <- tibble(pair = getPairIds(lnkmat), lnk = vectorize(lnkmat))
        
        
        cat('Saving file as:', paste0(output_dir, "/", cell, "_",d, "_xCellxSubj_reprod.rds"), '\n')
        saveRDS(lnktbl, file = paste0(output_dir, "/", cell, "_",d, "_xCellxSubj_reprod.rds"))
      })
    }
    
    
  }
  
  
}




# ============================================================================================================================================================
# ============================================================================================================================================================
# ============================================================================================================================================================
# FISHER TEST  
# ============================================================================================================================================================
# ============================================================================================================================================================
# ============================================================================================================================================================



compute_fisher_XcellXsubj <- function() {
  # ============================================================================
  # Adapted from Pan 
  # ============================================================================
  fisher <- function(predicted, notPredicted, trueSet, alternative) {
    # return 1. a fisher's test object, 2. contingency table
    
    predicted_false <- notPredicted
    predicted_true <- predicted %>% setdiff(predicted_false) # make sure the two are mutually exclusive
    background <- c(predicted_true, predicted_false)
    truth <- trueSet %>% intersect(background)
    
    a <- predicted_true %>% intersect(truth) 
    b <- predicted_true %>% setdiff(a)
    
    c <- predicted_false %>% intersect(truth)
    d <- predicted_false %>% setdiff(c)
    
    contingency <- data.frame(
      true = c(length(a), length(c)), 
      false = c(length(b), length(d)), 
      row.names = c("sampled", "not_sampled")
    )
    
    fisher <- fisher.test(contingency, alternative = alternative)
    
    # also return some basic stats
    # 1. or - odds ratio, 2. enrich, 3. depriv, 4. expected, ... etc
    
    contingency_null <- chisq.test(contingency)$expected
    
    stats <- c()
    stats["or"] <- (contingency[1, 1] / contingency[1, 2]) / (contingency[2, 1] / contingency[2, 2])
    stats["enrich"] <- (contingency[1, 1] / contingency[1, 2])
    stats["depriv"] <- (contingency[2, 1] / contingency[2, 2])
    stats["times_over_random"] <- (contingency[1, 1] / contingency_null[1, 1])
    stats["recovered_n"] <- contingency[1, 1]
    stats["recovered_frac"] <- contingency[1, 1] / length(trueSet) 
    
    return(list(tbl = contingency, tbl_null = contingency_null, test = fisher, stats = stats))
  }
  
  
  xsubj_fishers_test <- function() {
    cat("Running Fisher's exact tests for reproducibility...\n")
    
    library(tidyverse)
    library(parallel)
    library(data.table)
    dataset <- c('Ling-2024','ROSMAP-Mathys-2024','MSBB-2024', 'RADC-2024','HBCC-2024','SZBDMulti-2024',
                 'MSSM-2024','Rexach-2024','ROSMAP-erosion-2023','Kamath-NPH-2023',
                 'Hoffman-2023','ROSMAP-Mathys-2019','ROSMAP-Mathys-2023','ROSMAP-TREM2-2020',
                 'ROSMAP-DeJager-2024','ROSMAP-DeJager-2023','SEA-AD-MTG-2024','SEA-AD-DLPFC-2024')
    cell_types <- c('Mic', 'Oli', 'Opc', 'Exc', 'Inh', 'Ast')
    
    results <- list()
    
    
    # Define input/output directory
    output_dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher"
    input_dir <- "/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/lnkTables"
    
    
    reprod_all_files <- list.files(input_dir, pattern = "_xCellxSubj_reprod.rds", full.names = TRUE)
    
    for (cell in cell_types) {
      cat('Processing cell type:', cell, '\n')
      reprod_files <- reprod_all_files[grepl(cell, basename(reprod_all_files))]
      
      if (length(reprod_files)==0) { next }
      
      lnks_T <- list()  
      lnks_F <- list() 
      
      for (file in reprod_files) {
        cat('Starting:', file, '\n')
        tbl <- readRDS(file)  %>% as.data.table()
        dataset_name <- gsub("_xCellxSubj_reprod.rds", "", basename(file))
        cat('Processing:', dataset_name, '\n')
        lnks_T[[dataset_name]] <- filter(tbl, lnk == 1)$pair
        lnks_F[[dataset_name]] <- filter(tbl, lnk == 0)$pair
        rm(tbl)  
        gc()  
      }
      
      
      # # Separate linked (lnk == 1) and non-linked (lnk == 0) pairs
      # lnks_T <- lapply(lnktbls, function(tbl) filter(tbl, lnk == 1)$pair)
      # lnks_F <- lapply(lnktbls, function(tbl) filter(tbl, lnk == 0)$pair)
      
      
      fisher_test_wrapper <- function(predicted, notPredicted, trueSet) {
        fisher_results <- fisher(predicted, notPredicted, trueSet, alternative = "greater")
        return(fisher_results)
      }
      
      cat("Running Fisher's tests in parallel...\n")
      res <- setNames(mclapply(names(lnks_T), function(currDat_T) {
        datasets_S <- setdiff(names(lnks_T), currDat_T)
        
        setNames(lapply(datasets_S, function(currDat_S) {
          fisher(
            predicted = lnks_T[[currDat_S]], 
            notPredicted = lnks_F[[currDat_S]], 
            trueSet = lnks_T[[currDat_T]], 
            alternative = "greater"
          )
        }), datasets_S)
      }, mc.cores = 4), names(lnks_T))
      
      results[[cell]] <- res
    }
    saveRDS(results, file = file.path(output_dir, "xCellxSubj_all_cell_types_repro_fisher.rds"))
    cat("Fisher's test results saved!\n")
    
  }
  
  # Run the function
  xsubj_fishers_test()
}




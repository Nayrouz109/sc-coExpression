

#### this is to plot the aggregated cor for a given gene pair per dataset per cell type 

source("/home/nelazzabi/rewiring/scripts/functions/topCoExpr/visalize.R")
allGenes_xpr = fread("/home/nelazzabi/rewiring/data/coExpr/0all-human/0.expr/allGenes_exp_summary_logCPM.csv") %>% setDT()
allGenes_avg <- allGenes_xpr[, .(target_avg_expr_log_cpm = mean(avg_expr_log_cpm, na.rm = TRUE)), 
                             by = .(gene, cell_type)]

AstCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Ast_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
ExcCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Exc_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
InhCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Inh_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
MicCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Mic_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
OliCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Oli_tf_allTrg_avgRnkCoexpr_binding_lit.rds")
OpcCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/Opc_tf_allTrg_avgRnkCoexpr_binding_lit.rds")

geneA_interest = "MEF2C"
geneB_interest = "MAML2"
AstCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/perDataset/Ast_tf_allTrg_avg_perStudy_RnkCoexpr_binding_lit.rds") %>% filter((geneA == geneA_interest & geneB == geneB_interest) |(geneA == geneB_interest & geneB == geneA_interest)) 
MicCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/perDataset/Mic_tf_allTrg_avg_perStudy_RnkCoexpr_binding_lit.rds") %>% filter((geneA == geneA_interest & geneB == geneB_interest) |(geneA == geneB_interest & geneB == geneA_interest))
OliCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/perDataset/Oli_tf_allTrg_avg_perStudy_RnkCoexpr_binding_lit.rds") %>% filter((geneA == geneA_interest & geneB == geneB_interest) |(geneA == geneB_interest & geneB == geneA_interest))
ExcCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/perDataset/Exc_tf_allTrg_avg_perStudy_RnkCoexpr_binding_lit.rds") %>% filter((geneA == geneA_interest & geneB == geneB_interest) |(geneA == geneB_interest & geneB == geneA_interest))
InhCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/perDataset/Inh_tf_allTrg_avg_perStudy_RnkCoexpr_binding_lit.rds") %>% filter((geneA == geneA_interest & geneB == geneB_interest) |(geneA == geneB_interest & geneB == geneA_interest))
OpcCor = readRDS("/home/nelazzabi/rewiring/data/coExpr/0all-human/1.2TFallTargets-coExpr-perCell/perDataset/Opc_tf_allTrg_avg_perStudy_RnkCoexpr_binding_lit.rds") %>% filter((geneA == geneA_interest & geneB == geneB_interest) |(geneA == geneB_interest & geneB == geneA_interest))

cor_list <- list(
  Ast = AstCor,
  Exc = ExcCor,
  Inh = InhCor,
  Mic = MicCor,
  Oli = OliCor,
  Opc = OpcCor 
) 

rm(AstCor)
rm(MicCor)
rm(OliCor)
rm(ExcCor)
rm(InhCor)
rm(OpcCor)

#### visualization function ================ below =======================================================================================

plot_gene_corr_vs_expr <- function(expr_df, cor_list, geneA_interest, geneB_interest = NULL, size_col = NULL) {
  
  # 1. Subset the correlation dfs for geneA of interest only
  cor_df <- purrr::map2_df(cor_list, names(cor_list), ~ 
                             .x %>% 
                             filter(geneA == geneA_interest) %>% 
                             mutate(cell_type = .y))
  
  # 2. Add expression for geneBs
  cor_df <- cor_df %>%
    left_join(expr_df, by = c("geneB" = "gene", "cell_type" = "cell_type"))
  
  # 3. Flag highlight genes and assign size
  cor_df <- cor_df %>%
    mutate(
      highlight = ifelse(!is.null(geneB_interest) & geneB %in% geneB_interest, "firebrick", "black"),
      pt_alpha  = ifelse(highlight == "black", 0.6, 1),
      pt_size   = ifelse(highlight == "firebrick" & !is.null(size_col),
                         .data[[size_col]], 
                         NA_real_)
    )
  
  # 4. Plot
  p <- ggplot(cor_df, aes(x = target_avg_expr_log_cpm, y = avg_rank_norm_corr)) +
    
    # black dots (fixed size)
    geom_point(data = cor_df %>% filter(highlight == "black"),
               aes(color = highlight, alpha = pt_alpha), size = 3) +
    
    # red dots: scaled by size_col if given, else fixed bigger size
    geom_point(data = cor_df %>% filter(highlight == "firebrick"),
               aes(color = highlight, alpha = pt_alpha,
                   size = if (!is.null(size_col)) pt_size else NULL),
               inherit.aes = TRUE,
               show.legend = FALSE) +
    
    scale_color_identity() +
    scale_alpha_identity() +
    scale_size_continuous(range = c(4, 3)) +  # adjust scaling range for size_col
    facet_wrap(~cell_type, scales = "fixed") +
    labs(
      x = "Avg Expression (logCPM)",
      y = paste("Avg Rank Norm Corr with", geneA_interest),
      title = paste("Expression vs Correlation for", geneA_interest)
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1),
      axis.text   = element_text(size = 15, face = "bold"),
      axis.title  = element_text(size = 16),
      legend.position = "none"
    )
  
  return(p)
}

#### 1. visualization avg rank vs avg expression ================  =======================================================================================
### all gene pairs with TF are in black 
### "target" of interest in red 
plot_gene_corr_vs_expr(
  expr_df = allGenes_avg,
  cor_list = cor_list,
  geneA_interest = "MEF2C",
  geneB_interest = c("RUNX1"),   # highlight these
  size_col = "binding_score"           # scale size if available
)




#### 2. visualization cor ================ TF & Target =======================================================================================
### avg across all datasets 
plot_pairwise_cor <- function(cor_list, geneA_interest, geneB_interest) { 
  library(dplyr) 
  library(ggplot2) 
  library(purrr) 
  # Step 1: extract correlations for geneA-geneB pair across cell types 
  df <- map_dfr(names(cor_list), function(ct) { 
    cor_df <- cor_list[[ct]] 
  # filter for both possible orientations of the pair 
    cor_pair <- cor_df %>% 
    filter((geneA == geneA_interest & geneB == geneB_interest) | 
             (geneA == geneB_interest & geneB == geneA_interest) ) %>% 
    mutate(CellType = ct) 
    return(cor_pair) }
  ) 
  # Step 2: average correlation per cell type 
  avg_df <- df %>% group_by(CellType) %>% summarise(avg_rank_norm_corr = mean(avg_rank_norm_corr, na.rm = TRUE), .groups = "drop") 
  # Step 3: plot 
  p <- ggplot(avg_df, aes(x = CellType, y = avg_rank_norm_corr, group = 1)) + 
    geom_line(color = "firebrick", size = 1) + 
    geom_point(color = "firebrick", size = 3) + 
    theme_classic() + 
    theme( axis.text.x = element_text(angle = 0, hjust = 1), 
           axis.text = element_text(size = 15, face = "bold"), 
           axis.title = element_text(size = 16) ) + 
    labs( title = paste0("Avg correlation of ", geneA_interest, " & ", geneB_interest), x = "Cell type", y = "Average correlation" ) 
  return(list(data = avg_df, plot = p)) 
  }

### per dataset   
plot_pairwise_cor <- function(cor_list, geneA_interest, geneB_interest) {
  library(dplyr)
  library(ggplot2)
  library(purrr)
  
  # Step 1: extract correlations for geneA-geneB pair across cell types
  df <- map_dfr(names(cor_list), function(ct) {
    cor_df <- cor_list[[ct]]
    
    # filter for both orientations of the pair
    cor_pair <- cor_df %>%
      filter(
        (geneA == geneA_interest & geneB == geneB_interest) |
          (geneA == geneB_interest & geneB == geneA_interest)
      ) %>%
      mutate(CellType = ct)
    
    return(cor_pair)
  })
  
  df <- df %>%
    mutate(CellType = factor(CellType, levels = unique(CellType))) %>% 
    mutate(dataset_id = str_remove(dataset_id, "^[A-Za-z]+_"))
  
  p <- ggplot() +
    # dataset-specific points and lines (grey)
    geom_line(
      data = df,
      aes(x = CellType, y = rank_norm_corr, group = dataset_id
          ),
      color = "grey70",
      alpha = 0.7, 
      size = 0.8
    ) +
    geom_point(
      data = df,
      aes(x = CellType, y = rank_norm_corr),
      color = "grey40",
      alpha = 0.9, 
      size = 2 
    ) +
    # overlay average per cell type (one line)
    geom_line(
      data = df %>% distinct(CellType, .keep_all = TRUE),
      aes(x = CellType, y = avg_rank_norm_corr, group = 1),
      color = "firebrick",
      size = 1.5
    ) +
    geom_point(
      data = df %>% distinct(CellType, .keep_all = TRUE),
      aes(x = CellType, y = avg_rank_norm_corr),
      color = "firebrick",
      size = 4
    ) +
    labs(
      x = "Cell Type",
      y = "Rank Normalized Correlation",
      title = paste0(geneA_interest, "-", geneB_interest, " correlations")
    ) +
    theme_minimal(base_size = 14) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
  
  return(p)
}


plot_pairwise_cor(cor_list, "MEF2C", "MAML2") 
plot_pairwise_cor(cor_list, "MEF2C", "JAK2") 

plot_pairwise_cor(cor_list, "MEF2C", "MAML2")$plot

plot_pairwise_cor(cor_list, "NFIB", "ARHGAP6")$plot #*exc
plot_pairwise_cor(cor_list, "NFIB", "TACR1")$plot   #*inh
plot_pairwise_cor(cor_list, "NFIB", "CA13")$plot    #inh

plot_pairwise_cor(cor_list, "CUX1", "GAB1")$plot    #inh and exc 
plot_pairwise_cor(cor_list, "NFIA", "NGEF")$plot    #inh and exc 
















#### genes =============================================================================
### python script ********************
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import json

### 5% genes summary 
genes = "/home/nelazzabi/rewiring/data/genes/genes-listsV2.json"
with open(genes, "r") as fp:
  genes_set = json.load(fp)

flttype_dict = genes_set["less_stringent"]
all_genes = sorted(set().union(*flttype_dict.values())) #union of all genes 

# Create binary matrix: rows = cell types, columns = genes
binary_matrix = pd.DataFrame(
  0, index=flttype_dict.keys(), columns=all_genes
)

# Fill in 1s for presence
for celltype, genes in flttype_dict.items():
  binary_matrix.loc[celltype, genes] = 1

binary_matrix.T.to_csv('/home/nelazzabi/rewiring/data/genes/lessStringent_binary.csv', index=False)

### R script ********************
genes_binary <- read.csv("/home/nelazzabi/rewiring/data/genes/lessStringent_binary.csv")
genes_binary[, 1:6]

pheatmap::pheatmap(
  genes_binary,
  color = c("black", "blue"),  # Color for 0s and 1s
  cluster_rows = TRUE,  # Do not cluster rows
  cluster_cols = TRUE,  # Do not cluster columns
  legend = TRUE,  # Show legend
  main = "Gene Presence in Cell Types"
)

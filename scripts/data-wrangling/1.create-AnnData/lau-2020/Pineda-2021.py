#!/~/anaconda3/envs/homl3/bin/python

import os
import sys
import json
import anndata as ad
import pandas as pd

# ======================
# Paths
# ======================
input_dir        = "/cosmos/data/project-data/NW-rewiring/data/pans/0.prcsd-raw/h5ad"
output_dir       = "/cosmos/data/project-data/NW-rewiring/data/pans/0.prcsd-raw"
select_genes_cpm = "/cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/cpm"
select_genes_raw = "/cosmos/data/project-data/NW-rewiring/data/pans/1.slct-genes/raw"

os.makedirs(output_dir, exist_ok=True)
os.makedirs(select_genes_cpm, exist_ok=True)
os.makedirs(select_genes_raw, exist_ok=True)

# ======================
# Step 1: Load h5ad
# ======================
file_path = os.path.join(input_dir, "pineda.h5ad")
adata = ad.read_h5ad(file_path)


# ======================
# Step 2: Wrangle obs columns
# ======================
# Rename "sample" → "patientID"
adata.obs = adata.obs.rename(columns={"Sample_ID": "patientID"})

disease_map = {"ALS": "ALS", "FTLD": "FTLD",  "PN": "CTL"}
adata.obs["disease"] = adata.obs["Condition"].map(disease_map).astype("category")
adata.obs["diagnosis"] = "ALS_FTLD"

celltype_map = {
    # Astrocytes
    "Astrocyte": "Ast",
    
    # Endothelial
    "Endothelial": "Endo",
    "Mural": "Endo",  # optional, depends if you want to include pericytes
    
    # Excitatory neurons
    "Ex.L2_L3": "Exc",
    "Ex.L3_L5": "Exc",
    "Ex.L5": "Exc",
    "Ex.L5_L6": "Exc",
    "Ex.L5b_UMN_CT": "Exc",
    "Ex.L5b_UMN_PT": "Exc",
    "Ex.L6b": "Exc",
    
    # Inhibitory neurons
    "In.5HT3aR_VIPneg": "Inh",
    "In.5HT3aR_VIPpos": "Inh",
    "In.Basket_PVALB": "Inh",
    "In.Chandelier_PVALB": "Inh",
    "In.Rosehip": "Inh",
    "In.SST_NPYneg": "Inh",
    "In.SST_NPYpos": "Inh",
    
    # Microglia
    "Microglia": "Mic",
    
    # OPCs → Oligodendrocyte lineage
    "OPC": "Opc",
    "Oligodendrocyte": "Oli",
    
    # Other cell types you may want to ignore
    "Fibroblast": None,
    "T_Cell": None,
}

# Apply mapping
adata.obs["Cell_type"] = adata.obs["CellType"].map(celltype_map)



# ======================
# Step 3: Save wrangled file
# ======================
print("handling genes")
adata.var_names = adata.var['feature_name'].astype(str)
adata.var.index.name = None
adata.var_names_make_unique()

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/cleaning')
from cleanCellsnGenes import clean_ctmat

cleaned_adata = clean_ctmat(adata, gene_thr=0.02, sample_thr=0.02)

print("done with genes")
wrangled_path = os.path.join(output_dir, "Pineda-2021.h5ad")
cleaned_adata.write(wrangled_path)
print(f"Saved wrangled AnnData → {wrangled_path}")



# ======================
# Step 4: Gene filtering
# ======================
genes_sets = "/home/nelazzabi/rewiring/data/genes/genes-listsV2.json"
with open(genes_sets, "r") as fp:
    genes_set = json.load(fp)

values = list(genes_set["less_stringent"].values())
intersection_set = set(values[0])
for value in values[1:]:
    intersection_set &= set(value)
intersection_set = list(intersection_set)

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/filter-cpm')
from slcGenes_cpm import slcGenes_and_cpm_adata

# Raw (no CPM normalization)
filtered_raw = slcGenes_and_cpm_adata(
    wrangled_path,
    intersection_set,
    cpm_normalize=False,
    log_transform=False
)

# CPM-normalized
filtered_cpm = slcGenes_and_cpm_adata(
    wrangled_path,
    intersection_set,
    cpm_normalize=True,
    log_transform=False
)

# ======================
# Step 5: Split by canonical cell types
# ======================
valid_cell_types = ["Mic", "Opc", "Oli", "Inh", "Exc", "Ast"]

for subset_name, dataset, outdir in [
    ("raw", filtered_raw, select_genes_raw),
    ("cpm", filtered_cpm, select_genes_cpm)
]:
    for cell_type in valid_cell_types:
        adata_sub = dataset[dataset.obs["Cell_type"] == cell_type].copy()
        if adata_sub.n_obs == 0:
            continue
        outfile = os.path.join(outdir, f"{cell_type}_Pineda-2021.h5ad")
        adata_sub.write(outfile)
        print(f"Saved {outfile} ({adata_sub.n_obs} cells, {adata_sub.n_vars} genes)")


# ======================
# Step 3b: Summarize counts
# ======================
def summarize_counts(adata_obj, label):
    return (
        adata_obj.obs.groupby(["patientID", "Cell_type"])
        .size()
        .reset_index(name=label)
    )

# Summaries
summary_wrangled = summarize_counts(adata, "wrangled")
summary_raw = summarize_counts(filtered_raw, "raw_filtered")
summary_cpm = summarize_counts(filtered_cpm, "cpm_filtered")

# Merge all summaries into one table
summary_all = (
    summary_wrangled
    .merge(summary_raw, on=["patientID","Cell_type"], how="outer")
    .merge(summary_cpm, on=["patientID","Cell_type"], how="outer")
)

# fill only numeric cols
num_cols = ["wrangled","raw_filtered","cpm_filtered"]
summary_all[num_cols] = summary_all[num_cols].fillna(0).astype(int)

# Save to CSV
summary_out = os.path.join(output_dir, "Pineda-2021_cell_counts_all2.csv")
summary_all.to_csv(summary_out, index=False)
print(f"Saved combined cell counts summary → {summary_out}")





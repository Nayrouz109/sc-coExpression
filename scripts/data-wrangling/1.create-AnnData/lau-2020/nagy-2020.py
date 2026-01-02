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
print("reading raw file")
file_path = os.path.join(input_dir, "nagy.h5ad")
adata = ad.read_h5ad(file_path)
adata.obs = adata.obs.reset_index().rename(columns={"index": "cell_id"}) #assign a cell ID column
print("done reading raw file")

# ======================
# Step 2: Wrangle obs columns
# ======================
# Rename "sample" → "patientID"
adata.obs["patientID"] = (
    adata.obs["cell_id"].str.extract(r"\.(\d+)").iloc[:, 0] + "_" +
    adata.obs["batch"].astype(str) + "_" +
    adata.obs["disease"].astype(str)
)

# Remap cell types
# Mapping dictionary
celltype_map = {
    # Astrocytes
    "Astros_1": "Ast",
    "Astros_2": "Ast",
    "Astros_3": "Ast",
    
    # Endothelial
    "Endo": "Endo",
    
    # Excitatory neurons
    "Ex_1_L5_6": "Exc",
    "Ex_2_L5": "Exc",
    "Ex_3_L4_5": "Exc",
    "Ex_4_L_6": "Exc",
    "Ex_5_L5": "Exc",
    "Ex_6_L4_6": "Exc",
    "Ex_7_L4_6": "Exc",
    "Ex_8_L5_6": "Exc",
    "Ex_9_L5_6": "Exc",
    "Ex_10_L2_4": "Exc",
    
    # Inhibitory neurons
    "Inhib_1": "Inh",
    "Inhib_2_VIP": "Inh",
    "Inhib_3_SST": "Inh",
    "Inhib_4_SST": "Inh",
    "Inhib_5": "Inh",
    "Inhib_6_SST": "Inh",
    "Inhib_7_PVALB": "Inh",
    "Inhib_8_PVALB": "Inh",
    
    # Microglia / Macrophages
    "Micro/Macro": "Mic",
    
    # OPCs → Oligodendrocyte precursors
    "OPCs_1": "Opc",
    "OPCs_2": "Opc",
    
    # Mature oligodendrocytes
    "Oligos_1": "Oli",
    "Oligos_2": "Oli",
    "Oligos_3": "Oli",
    
    # Any "Mix" clusters could be ignored or mapped to NaN if you don't want them
    "Mix_1": None,
    "Mix_2": None,
    "Mix_3": None,
    "Mix_4": None,
    "Mix_5": None,
}

adata.obs["Cell_type"] = adata.obs["cell_type"].map(celltype_map)
adata.obs["Cell_type"] = adata.obs["Cell_type"].astype("category")
print("done mapping cell types")

disease_map = {"MDD": "MDD", "Control": "CTL"}
adata.obs["disease"] = adata.obs["disease"].map(disease_map).astype("category")
adata.obs["diagnosis"] = "MDD"
print("done mapping diagnosis")

# ======================
# Step 3: Save wrangled file
# ======================
print("changing gene names")
adata.var_names = adata.var['hgnc']
adata.var['ensembl_id'] = adata.var['feature_id']
print("done changing gene names")

sys.path.append('/home/nelazzabi/rewiring/scripts/functions/cleaning')
from cleanCellsnGenes import clean_ctmat

cleaned_adata = clean_ctmat(adata, gene_thr=0.02, sample_thr=0.02)

print("saving")
wrangled_path = os.path.join(output_dir, "Nagy-2020.h5ad")
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
        outfile = os.path.join(outdir, f"{cell_type}_Nagy-2020.h5ad")
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
summary_out = os.path.join(output_dir, "Nagy-2020_cell_counts_all2.csv")
summary_all.to_csv(summary_out, index=False)
print(f"Saved combined cell counts summary → {summary_out}")







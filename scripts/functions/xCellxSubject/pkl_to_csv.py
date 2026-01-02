




###============================================================
### 1. reading in the pkl in Python 
###============================================================
import pickle
import pandas as pd
import os


with open("/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher/xCellxSubj_all_cell_types_repro_fisher.pkl", "rb") as f:
    results = pickle.load(f)

# Example: list all cell types
print(results.keys())

# Example: look at one cell type
mic_res = results["Mic"]

# Example: look at comparisons for DatasetA
print(mic_res["DatasetA"].keys())

# Example: get Fisher result of DatasetB vs DatasetA
res_AB = mic_res["DatasetA"]["DatasetB"]
print(res_AB)



###============================================================
### 1. reading in the pkl in Python 
###============================================================


# Flatten results into a long-format dataframe
flat_records = []
for cell_type, datasets in results.items():
    for reference_dataset, preds in datasets.items():
        for dataset, fisher_res in preds.items():
            stats = fisher_res["stats"]
            flat_records.append({
                "reference_dataset": reference_dataset.replace(f"{cell_type}_", ""),
                "dataset": dataset.replace(f"{cell_type}_", ""),
                "cell_type": cell_type,
                "odds_ratio": stats.get("or", None),
                "pvalue": stats.get("p_value", None),
                "reprod": "xcellxSubj"
            })

df = pd.DataFrame(flat_records)

# Save to CSV
output_dir = "/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/fisher"
csv_file = os.path.join(output_dir, "xCellxSubj_all_cell_types_repro_fisher.csv")
df.to_csv(csv_file, index=False)

print(f"✅ Flattened results saved to {csv_file}")
print(df.head())

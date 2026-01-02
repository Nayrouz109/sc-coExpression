

import os
import glob
import pandas as pd

# Path to your directory
path = "/cosmos/data/project-data/NW-rewiring/coExpr/xCellxSubject/rpr/lnkTables"

# Get only the parquet files that start with "Mic"
files = sorted(glob.glob(os.path.join(path, "Mic*_xcellxSubject_reprod.parquet")))



# Read all files into dataframes
dfs = {os.path.basename(f): pd.read_parquet(f, columns=["pair"]) for f in files[1:4]}

# 1. Check that the number of unique values in pair_ids are the same
unique_counts = {fname: df["pair"].nunique() for fname, df in dfs.items()}
print("\nUnique counts of pair_ids:")
for k, v in unique_counts.items():
    print(f"{k}: {v}")

if len(set(unique_counts.values())) == 1:
    print("\n✅ All files have the same number of unique pair_ids.")
else:
    print("\n⚠️ Files differ in the number of unique pair_ids!")

# 2. Check that the order of pair_ids is the same
pair_id_orders = [df["pair"].tolist() for df in dfs.values()]

# Compare each file's order with the first file
first_file = list(dfs.keys())[0]
first_order = pair_id_orders[0]

order_mismatches = {}
for fname, order in zip(list(dfs.keys()), pair_id_orders):
    if order != first_order:
        order_mismatches[fname] = True

if not order_mismatches:
    print("\n✅ The order of pair_ids is the same across all files.")
else:
    print("\n⚠️ The order of pair_ids differs in these files:")
    for fname in order_mismatches:
        print("-", fname)

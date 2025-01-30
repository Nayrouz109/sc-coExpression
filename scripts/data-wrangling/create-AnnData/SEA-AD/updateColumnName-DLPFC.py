

#!~/anaconda3/envs/homl3/bin/python

import anndata as ad
import os
import glob


data_dir   = '/cosmos/data/downloaded-data/AD-meta/SEA-AD-DLPFC/data/0-source/expr/annData-cellTpe/'
files = glob.glob(os.path.join(data_dir, "*.h5ad"))

# Iterate over each file, load it, update the column name, and save it back
for file_path in files:
    file_name = os.path.basename(file_path)  
    print(f"Processing file: {file_name}")
    
    data = ad.read_h5ad(file_path)
    print(f"Loaded {file_name}")

    # Check if "Donor ID" exists in data.obs before renaming
    if "Donor ID" in data.obs.columns:
        # Rename "Donor ID" column to "patientID"
        data.obs = data.obs.rename(columns={"Donor ID": "patientID"})
        print(f"Renamed 'Donor ID' to 'patientID' in {file_name}")
    else:
        print(f"'Donor ID' column not found in {file_name}")

    # Save the updated AnnData object, overwriting the original file
    data.write_h5ad(file_path, compression="gzip")
    print(f"Saved updated file: {file_name}\n")

print("All files processed successfully.")

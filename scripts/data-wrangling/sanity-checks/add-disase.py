#!~/anaconda3/envs/homl3/bin/python

import os
import pandas as pd
import anndata as ad
from pathlib import Path

# Extract study name from filename
def extract_metadata_from_filename(filename):
    """Extracts the study name from the filename."""
    name_parts = filename.split('_')
    cell_type = name_parts[0] 
    study_name = name_parts[-1].split('.')[0] 
    return cell_type, study_name


# Directories
input_dir = Path("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag")
output_dir = Path("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease")
output_dir.mkdir(parents=True, exist_ok=True)

# Load the big metadata file
p_meta = pd.read_csv("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData-updated.csv")
unique_studies = p_meta["study"].unique()

# Initialize results
missing_study_files = []
missing_patient_files = []
summary_data = []

# Loop over all .h5ad files in the input directory
for file_path in input_dir.glob("*.h5ad"):
    print(f"\nProcessing file: {file_path.name}")
    
    try:
        # Load AnnData object
        adata = ad.read_h5ad(file_path)
        adata.obs["patientID"] = adata.obs["patientID"].astype(str)
    except Exception as e:
        print(f"Error reading {file_path.name}: {e}")
        continue
    
    # Extract study info from annData  
    cell_type, study_name = extract_metadata_from_filename(file_path.name)
    ad_patient_ids = adata.obs["patientID"].astype(str).unique()
    
    # Extract study info from big patient meta file 
    p_meta_study = p_meta[p_meta["study"] == study_name].drop_duplicates(subset=["patientID"])
    p_meta_study_patient_ids = p_meta_study["patientID"].astype(str).unique()

    # Check if study name exists in the metadata
    if study_name not in unique_studies:
        missing_study_files.append(file_path.name)
        print(f"Study name {study_name} is missing in the patient meta file.")
        continue
    
    # Check for missing patient IDs
    missing_patient_ids = [pid for pid in ad_patient_ids if pid not in p_meta_study_patient_ids]
    if missing_patient_ids:
        missing_patient_files.append(file_path.name)
        print(f"Missing patient IDs in {file_path.stem}: {missing_patient_ids}")
    
    # Integrate diagnosis and disease information
    adata.obs = adata.obs.merge(
        p_meta_study[["patientID", "disease"]], 
        on="patientID", 
        how="left") 
    
    # Ensure diagnosis values are valid
    valid_diagnoses = {"yes", "no", pd.NA}
    observed_diagnoses = set(adata.obs["diagnosis"].dropna().unique())
    if not valid_diagnoses.issuperset(observed_diagnoses) or {"yes", "no"}.difference(observed_diagnoses):
        raise ValueError(f"Invalid or insufficient diagnosis values in {file_path.name}: {observed_diagnoses}")
    
    # Summarize file information
    summary_data.append({
        "study name": study_name,
        "cell type": cell_type,
        "total genes": adata.shape[1],
        "total single cells": adata.shape[0],
        "total patient IDs": len(ad_patient_ids), 
        "missing patient IDs": len(missing_patient_ids)  # Include the count of missing IDs in the summary
    })
    
    # Save updated AnnData
    output_file_path = output_dir / file_path.name
    try:
        adata.write_h5ad(output_file_path, compression="gzip")
        print(f"Integrated diagnosis and disease into {file_path.stem} and saved to {output_file_path}")
    except Exception as e:
        print(f"Error saving {file_path.name}: {e}")

# Save the summary to a CSV file
summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(output_dir / "summary.csv", index=False)
print("\nSummary saved to summary.csv")

# Print results
print(f"\nFiles with missing study names: {missing_study_files}")
print(f"Files with missing patient IDs: {missing_patient_files}")


import os
import anndata as ad
import numpy as np
import pandas as pd
import random

def downsample_anndata(adata, n_samples):
    """Downsample the AnnData object to `n_samples` cells per sample."""
    if adata.n_obs > n_samples:
        sampled_indices = np.random.choice(adata.n_obs, n_samples, replace=False)
        return adata[sampled_indices, :]
    return adata


def process_files(input_dir, output_dir_samples, output_dir_cells, cell_type, conditions):
    os.makedirs(output_dir_samples, exist_ok=True)
    os.makedirs(output_dir_cells, exist_ok=True)
    
    # Identify relevant files
    files = [f for f in os.listdir(input_dir) if f.startswith(cell_type) and f.endswith(".h5ad")]
    
    for file in files:
        print(f"\nProcessing file: {file}\n")
        file_path = os.path.join(input_dir, file)
        adata = ad.read_h5ad(file_path)
        
        # 0.1 Ensure the dataset has 'patientID', and 'disease' columns in metadata
        if 'disease' not in adata.obs.columns or 'patientID' not in adata.obs.columns:
            print(f"Skipping {file}: Missing required metadata columns.")
            continue
        
        # 0.2 Check if disease categories exist in metadata
        available_categories = set(adata.obs["disease"].unique())
        if not set(conditions).issubset(available_categories):
            print(f"Skipping {file_path}: Required categories {conditions} not found ({available_categories})")
            continue
        
        # 0.3 Remove patientIDs with fewer than 20 single cells
        patient_cell_counts = adata.obs.groupby("patientID").size()
        valid_patients = patient_cell_counts[patient_cell_counts >= 20].index
        adata = adata[adata.obs["patientID"].isin(valid_patients)] 


        # Step 1: Downsample samples per condition if necessary
        print("downsampling patientIDs")
        grouped = adata.obs.groupby("disease")["patientID"].nunique()
        min_samples = grouped.min()
        
        if grouped.nunique() == 1:
            print(f"Conditions already balanced for {file}. Skipping sample downsampling.")
            adata_sampled = adata  # Keep as is
            output_file_samples = os.path.join(output_dir_samples, file)
            adata_sampled.write_h5ad(output_file_samples, compression="gzip")
        else:
            sampled_samples = []
            for cond in conditions:
                condition_samples = adata.obs[adata.obs["disease"] == cond]["patientID"].unique()
                sampled_samples.extend(random.sample(list(condition_samples), min_samples))
            adata_sampled = adata[adata.obs["patientID"].isin(sampled_samples)]
            
            # Save downsampled per sample
            output_file_samples = os.path.join(output_dir_samples, file)
            adata_sampled.write_h5ad(output_file_samples, compression="gzip")
        
        # Step 2: Downsample single cells per sample
        print("downsampling cells")
        cell_counts = adata_sampled.obs.groupby("patientID").size()
        min_cells = cell_counts.min()
        
        downsampled_cells = []
        for sample in adata_sampled.obs["patientID"].unique():
            sub_adata = adata_sampled[adata_sampled.obs["patientID"] == sample]
            downsampled_cells.append(downsample_anndata(sub_adata, min_cells))
        
        adata_final = ad.concat(downsampled_cells, merge="same")
        
        # Save final downsampled dataset
        output_file_cells = os.path.join(output_dir_cells, file)
        adata_final.write_h5ad(output_file_cells , compression= "gzip")
         
        print(f"Processed and saved: {file}")

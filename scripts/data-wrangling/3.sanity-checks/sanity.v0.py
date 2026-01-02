
#!~/anaconda3/envs/homl3/bin/python
import numpy as np
import pandas as pd
import anndata as ad
import os
from scipy.sparse import csr_matrix



# Define the directory containing the files
directory = '/space/scratch/nairuz-rewiring/data/0.prcsd-raw'

# Function to check if "patientID" column is present in the .obs attribute
def check_patient_id_in_obs(file_path):
    # Load the .h5ad file using Anndata
    try:
        print(f"checking file: {filename}")
        adata = ad.read_h5ad(file_path)
        # Check if 'patientID' column exists in .obs
        if 'patientID' not in adata.obs.columns:
            return True  # 'patientID' is missing
        return False  # 'patientID' exists
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return False

# List to store filenames lacking 'patientID' in .obs
missing_patient_id_files = []

# List all files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.h5ad'):
        file_path = os.path.join(directory, filename)
        if check_patient_id_in_obs(file_path):
            missing_patient_id_files.append(filename)

# Print the list of files lacking 'patientID'
if missing_patient_id_files:
    print("The following files lack the 'patientID' column in .obs:")
    for file in missing_patient_id_files:
        print(file)
else:
    print("All .h5ad files have the 'patientID' column in .obs.")

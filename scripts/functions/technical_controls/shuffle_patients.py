

import numpy as np
import scanpy as sc  # Assuming AnnData is being used

def filter_and_shuffle_patients(adata, categories):
    """
    Filters patient IDs based on disease categories, shuffles them into two random groups,
    and assigns new disease labels.
    
    Parameters:
    adata (AnnData): Annotated data matrix.
    categories (list): List of disease categories to filter patients.
    seed (int): Random seed for reproducibility.
    
    Returns:
    AnnData: A new AnnData object with updated disease labels.
    """
    # Filter patient IDs based on disease categories
    valid_patients = adata.obs[adata.obs["disease"].isin(categories)]["patientID"].unique()
    filtered_data = adata[adata.obs["patientID"].isin(valid_patients)].copy()
    
    # Shuffle patients into two random groups
    all_patients = list(valid_patients)
    np.random.shuffle(all_patients)
    midpoint = len(all_patients) // 2
    grp1_patients = set(all_patients[:midpoint])
    grp2_patients = set(all_patients[midpoint:])
    
    # Assign random group labels
    filtered_data.obs["random_group"] = filtered_data.obs["patientID"].map(lambda x: "grp1" if x in grp1_patients else "grp2")
    filtered_data.obs = filtered_data.obs.rename(columns={
        "disease": "old_disease",
        "random_group": "disease"
    })
    
    # Replace group labels with new disease labels
    filtered_data.obs["disease"] = filtered_data.obs["disease"].replace({"grp1": "AD", "grp2": "CTL"})
    
    return filtered_data

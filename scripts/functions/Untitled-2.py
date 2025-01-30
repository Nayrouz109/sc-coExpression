

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import bottleneck 
from sklearn.preprocessing import QuantileTransformer
from scipy.sparse import issparse


def rank(data):
    """Rank normalize data
    
    Rank standardize data to make nonparametric
    
    Arguments:
        data {np.array} -- 2-D coexpression network
    
    Returns:
        np.array -- Rank normalized between 0 and 1 array
    """
    orig_shape = data.shape
    data = bottleneck.nanrankdata(data) - 1
    return (data / np.sum(~np.isnan(data))).reshape(orig_shape)

def compute_gene_correlation(adata):
    """
    Compute the Pearson correlation coefficients, rank normalize them,
    and aggregate by patientID and cell type.
    
    Arguments:
        adata {AnnData} -- The input AnnData object containing gene expression data and metadata
    
    Returns:
        dict -- Dictionary containing:
            - 'correlation_per_patientID' : raw Pearson correlation values per patientID
            - 'rank_normalized_per_patientID' : rank normalized correlation matrix per patientID
            - 'aggregated_rank_normalized_matrix' : rank normalized aggregated correlation matrix per cell type
    """
    # Step 1: Subset the data by patientID (based on metadata)
    patient_ids = adata.obs['patientID'].unique()
    
    # Initialize dictionaries to hold the results
    correlation_per_patientID = {}
    rank_normalized_per_patientID = {}
    aggregated_rank_normalized_matrix = {}

    # For each patient, compute the Pearson correlation matrix
    for patient_id in patient_ids:
        # Subset the expression matrix for the current patient
        subset_data = adata[adata.obs['patientID'] == patient_id].X

        if issparse(subset_data):
            subset_data = subset_data.toarray()
        
        if subset_data.shape[0] > 19:
            # Compute the Pearson correlation matrix
            corr_matrix = np.corrcoef(subset_data.T)
                
        # Store raw correlation values
        correlation_per_patientID[patient_id] = corr_matrix
        
        # Step 2: Rank normalize the correlation matrix
        rank_normalized_matrix = rank(corr_matrix)
        rank_normalized_per_patientID[patient_id] = rank_normalized_matrix
        
    # Step 3: Aggregate the rank-normalized matrices by cell type
    # Group by cell type and calculate the mean rank-normalized correlation
    cell_types = adata.obs['CellType'].unique()
    
    for cell_type in cell_types:
        # Get all rank-normalized matrices for the current cell type across patients
        cell_type_data = []
        for patient_id in patient_ids:
            patient_subset = adata[(adata.obs['CellType'] == cell_type) & (adata.obs['patientID'] == patient_id)].X
            corr_matrix = np.corrcoef(patient_subset.T)
            rank_normalized_matrix = rank(corr_matrix)
            cell_type_data.append(rank_normalized_matrix)
        
        # Average all rank-normalized matrices for the current cell type
        aggregated_matrix = np.mean(cell_type_data, axis=0)
        aggregated_rank_normalized_matrix[cell_type] = aggregated_matrix
    
    return {
        'correlation_per_patientID': correlation_per_patientID,
        'rank_normalized_per_patientID': rank_normalized_per_patientID,
        'aggregated_rank_normalized_matrix': aggregated_rank_normalized_matrix
    }








test = p_meta[p_meta["study"] == "ROSMAP-DeJager-2024" ]
test = p_meta[p_meta["study"] == "ROSMAP-Mathys-2019" ]

adata = ad.read_h5ad("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/Mic_ROSMAP-Mathys-2023.h5ad")
adata = ad.read_h5ad("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/Mic_ROSMAP-Mathys-2023.h5ad")

test_patient_ids = test["patientID"].unique()
adata_patient_ids = adata.obs["patientID"].unique()
pd.Series(adata_patient_ids).isin(test_patient_ids).all()
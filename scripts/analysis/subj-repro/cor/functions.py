
import numpy as np
import scipy.sparse
import pandas as pd 
import anndata as ad 
import os 
from anndata import read_h5ad
from joblib import Parallel, delayed



# Function to compute rank normalization
def rank(data):
    """Rank normalize data
    
    Rank standardize data to make nonparametric
    
    Arguments:
        data {np.array} -- 2-D coexpression network
    
    Returns:
        np.array -- Rank normalized between 0 and 1 array
    """
    orig_shape = data.shape
    #data = bottleneck.nanrankdata(data) - 1
    #return (data / np.sum(~np.isnan(data))).reshape(orig_shape)
    data = np.argsort(np.argsort(data, axis=None))  # Rank data
    return (data / (data.size - 1)).reshape(orig_shape)


# Function to compute sparse Pearson correlation
def sparse_corr(A):
    """Compute Pearson correlation for a sparse matrix.
    
    Arguments:
        A {scipy.sparse matrix} -- Sparse matrix of gene expression values
    Returns:
        np.matrix -- Correlation matrix
    """
    N = A.shape[0]
    #C = ((A.T * A - (A.sum(axis=0).T * A.sum(axis=0) / N)) / (N - 1)).todense()
    C = ((A.T*A -(sum(A).T*sum(A)/N))/(N-1)).todense()
    V = np.sqrt(np.mat(np.diag(C)).T * np.mat(np.diag(C)))
    COR = np.divide(C, V + 1e-119)
    return COR

# Function to process individual patient

def process_patient(patient_id, adata, min_cells_threshold, correlation_type, replace_nans):
    print(f"processing {patient_id}\n")
    subset_data = adata[adata.obs['patientID'] == patient_id].X
    subset_data = scipy.sparse.csr_matrix(subset_data)

    if subset_data.shape[0] < min_cells_threshold:
        print(f"Skipping patient {patient_id} due to insufficient number of cells.\n")
        return None, None

    # Compute correlation matrix
    if correlation_type == 'pearson':
        corr_matrix = sparse_corr(subset_data)
    else:
        raise ValueError(f"Unsupported correlation type: {correlation_type}")

    np.fill_diagonal(corr_matrix, 1)
    rank_normalized_matrix = rank(corr_matrix)

    if replace_nans:
        rank_normalized_matrix[np.isnan(rank_normalized_matrix)] = np.nanmean(rank_normalized_matrix)

    gene_names = adata.var_names.to_list()
    print(f"Finished processing patient {patient_id}\n")

    return corr_matrix, rank_normalized_matrix, gene_names

# Process datasets and save results
# Parameters
min_cells_threshold = 20
correlation_type = "pearson"
replace_nans = True

def process_and_save(file_path, patient_ids, study_name, output_dirs, brain_region=None):
    print(f"Loading data from {file_path}\n")
    adata = read_h5ad(file_path)

    # Subset data to specific patient IDs and brain region if applicable
    adata = adata[adata.obs['patientID'].isin(patient_ids)]
    if brain_region:
        adata = adata[adata.obs['region'] == brain_region]

    summary_data = []

    results = Parallel(n_jobs=7)(
        delayed(process_patient)(patient_id, adata, min_cells_threshold, correlation_type, replace_nans)
        for patient_id in patient_ids
    )

    for patient_id, (corr_matrix, rank_matrix, gene_names) in zip(patient_ids, results):
        if corr_matrix is None or rank_matrix is None or gene_names is None:
            continue

        # Define metadata
        #cell_type = "Exc"
        cell_type = "Opc"
        output_base = f"{cell_type}_{patient_id}_{study_name}_{brain_region or 'NA'}"

        # Save raw correlation matrix
        raw_path = os.path.join(output_dirs["raw"], f"{output_base}_raw.h5ad")
        raw_adata = ad.AnnData(X=corr_matrix)
        raw_adata.var_names = gene_names
        raw_adata.write_h5ad(raw_path, compression="gzip")

        # Save rank-normalized matrix
        rank_path = os.path.join(output_dirs["rank"], f"{output_base}_rank.h5ad")
        rank_adata = ad.AnnData(X=rank_matrix)
        rank_adata.var_names = gene_names
        rank_adata.write_h5ad(rank_path, compression="gzip")

        # Add to summary grid data
        summary_data.append({
            "PatientID": patient_id,
            "StudyName": study_name,
            "BrainRegion": brain_region or 'NA',
            "RawMatrixPath": raw_path,
            "RankMatrixPath": rank_path
        })

        print(f"Saved results for patient {patient_id}\n")

    # Convert summary data to DataFrame
    summary_df = pd.DataFrame(summary_data)

    # Save summary grid
    summary_path = os.path.join(output_dirs["raw"], f"summary_{study_name}_{brain_region or 'NA'}.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary grid to {summary_path}\n")

    return summary_df
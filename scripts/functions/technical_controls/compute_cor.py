

from pathlib import Path
import sys
import gc 


sys.path.append('/home/nelazzabi/rewiring/scripts/functions/compute-cor')
from corOnly_parallelV6_disease import process_patient


def compute_gene_correlation(adata, 
                             correlation_type='pearson', 
                             replace_nans=True, 
                             min_cells_threshold=20, 
                             category=None): 
    print("Starting correlation computation...")

    print(f"Processing subset for category: {category}")

    subset = adata[adata.obs['disease'] == category]
    patient_ids = subset.obs['patientID'].unique()
    total_patients = len(patient_ids)
    print(f"Total number of patients in {category}: {total_patients}")

    aggregated_matrix = None
    gene_names = None
    patient_count = 0

    for patient_id in patient_ids:
        print(f"Processing patient ID: {patient_id}")

        result = process_patient(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)

        if result is not None and result[1] is not None:
            patient_matrix = result[1]["rank_normalized_matrix"]
            if gene_names is None:
                gene_names = result[1]["gene_names"]

            # Aggregate the patient matrix into the aggregated matrix
            if aggregated_matrix is None:
                aggregated_matrix = patient_matrix
            else:
                aggregated_matrix += patient_matrix

            patient_count += 1
            print(f"Patient {patient_id} matrix aggregated.")
            del patient_matrix
            gc.collect() 
        else:
            print(f"No valid matrix for patient {patient_id}.")
    
    # Take the average after all patient matrices have been added
    if aggregated_matrix is not None and patient_count > 0:
        aggregated_matrix = aggregated_matrix / patient_count

    category_data =  aggregated_matrix 
    #category_data = {
        #'aggregated_rank_normalized_matrix': aggregated_matrix,
        #'gene_names': gene_names
    #}

    return(category_data)




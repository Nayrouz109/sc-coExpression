

import numpy as np
import pandas as pd
import anndata as ad 
from scipy.stats import fisher_exact

def binarize_network(adj_matrix, top_percentile=1):
    """
    Binarize the adjacency matrix by selecting the top X% edges as positive.

    Parameters:
        adj_matrix (numpy.ndarray): Adjacency matrix of the co-expression network.
        top_percentile (float): Top percentile to consider as positive edges.

    Returns:
        binarized_matrix (numpy.ndarray): Binary adjacency matrix with 1s for positive edges and 0s for negative.
    """
    threshold = np.percentile(adj_matrix, 100 - top_percentile)
    binarized_matrix = (adj_matrix >= threshold).astype(int)
    return binarized_matrix

def compute_odds_ratio(ref_matrix, test_matrix):
    """
    Compute the odds ratio and Fisher's exact test p-value between two networks.

    Parameters:
        ref_matrix (numpy.ndarray): Binarized adjacency matrix of the reference network.
        test_matrix (numpy.ndarray): Binarized adjacency matrix of the test network.

    Returns:
        odds_ratio (float): The odds ratio value.
        p_value (float): The p-value from Fisher's exact test.
    """
    # Convert adjacency matrices to boolean
    ref_positive = ref_matrix == 1
    test_positive = test_matrix == 1

    ref_negative = ref_matrix == 0
    test_negative = test_matrix == 0

    # Compute overlaps
    pos_pos_overlap = np.sum(ref_positive & test_positive)
    pos_neg_overlap = np.sum(ref_positive & test_negative)
    neg_pos_overlap = np.sum(ref_negative & test_positive)
    neg_neg_overlap = np.sum(ref_negative & test_negative)

    # Construct the contingency table
    contingency_table = np.array([
        [pos_pos_overlap, pos_neg_overlap],
        [neg_pos_overlap, neg_neg_overlap],
    ])

    # Compute odds ratio and Fisher's exact test p-value
    odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")

    return odds_ratio, p_value



def main(ref_adj_matrix, test_adj_matrix, top_percentile=1):
    """
    Main function to calculate odds ratio and p-value for reproducibility.

    Parameters:
        ref_adj_matrix (numpy.ndarray): Adjacency matrix of the reference network.
        test_adj_matrix (numpy.ndarray): Adjacency matrix of the test network.
        top_percentile (float): Top percentile to consider as positive edges.

    Returns:
        result (dict): A dictionary containing the odds ratio and p-value.
    """
    # Binarize the networks
    ref_binarized = binarize_network(ref_adj_matrix, top_percentile)
    test_binarized = binarize_network(test_adj_matrix, top_percentile)

    # Compute odds ratio and p-value
    odds_ratio, p_value = compute_odds_ratio(ref_binarized, test_binarized)

    return {
        "odds_ratio": odds_ratio,
        "p_value": p_value
    }


def compute_avg_odds_ratio(ref_matrix_path, test_matrix_path, top_percentile=1):
    # Load adjacency matrices
    ref_adj_matrix = ad.read_h5ad(ref_matrix_path).X
    test_adj_matrix = ad.read_h5ad(test_matrix_path).X

    # Symmetrize matrices
    ref_adj_matrix = np.triu(ref_adj_matrix, 1) + np.triu(ref_adj_matrix, 1).T
    test_adj_matrix = np.triu(test_adj_matrix, 1) + np.triu(test_adj_matrix, 1).T

    # Compute odds ratio twice with alternating reference and test designations
    odds_ratio_1 = main(ref_adj_matrix, test_adj_matrix, top_percentile)
    odds_ratio_2 = main(test_adj_matrix, ref_adj_matrix, top_percentile)

    # Return average odds ratio
    return np.mean([odds_ratio_1["odds_ratio"], odds_ratio_2["odds_ratio"]])



def perform_comparisons(grid):
    comparisons = []

    # Comparison 1: Same patient, different studies
    print("Starting Comparison 1: Same patient, different studies\n")
    for patient_id in grid['PatientID'].unique():
        patient_data = grid[grid['PatientID'] == patient_id]
        if len(patient_data['StudyName'].unique()) > 1:
            studies = patient_data['StudyName'].unique()
            for i, study1 in enumerate(studies):
                for study2 in studies[i + 1:]:
                    print(f"Comparing patient {patient_id} between studies {study1} and {study2}\n")
                    ref_path = patient_data[patient_data['StudyName'] == study1]['RankMatrixPath'].values[0]
                    test_path = patient_data[patient_data['StudyName'] == study2]['RankMatrixPath'].values[0]
                    avg_odds_ratio = compute_avg_odds_ratio(ref_path, test_path)
                    comparisons.append((patient_id, study1, study2, avg_odds_ratio, "Same patient, different studies"))

    print("Finished Comparison 1\n")

    # Comparison 2: Random patient, different studies
    print("Starting Comparison 2: Random patient, different studies\n")
    for i, row1 in grid.iterrows():
        for j, row2 in grid.iterrows():
            if i != j and row1['StudyName'] != row2['StudyName']:
                print(f"Comparing random patients {row1['PatientID']} and {row2['PatientID']} between studies {row1['StudyName']} and {row2['StudyName']}\n")
                avg_odds_ratio = compute_avg_odds_ratio(row1['RankMatrixPath'], row2['RankMatrixPath'])
                comparisons.append((
                    f"{row1['PatientID']} vs {row2['PatientID']}",
                    row1['StudyName'],
                    row2['StudyName'],
                    avg_odds_ratio,
                    "Random patient, different studies"
                ))

    print("Finished Comparison 2\n")

    # Comparison 3: Random patient, same study
    print("Starting Comparison 3: Random patient, same study\n")
    for i, row1 in grid.iterrows():
        for j, row2 in grid.iterrows():
            if i != j and row1['StudyName'] == row2['StudyName']:
                print(f"Comparing random patients {row1['PatientID']} and {row2['PatientID']} within study {row1['StudyName']}\n")
                avg_odds_ratio = compute_avg_odds_ratio(row1['RankMatrixPath'], row2['RankMatrixPath'])
                comparisons.append((
                    f"{row1['PatientID']} vs {row2['PatientID']}",
                    row1['StudyName'],
                    row2['StudyName'],
                    avg_odds_ratio,
                    "Random patient, same study"
                ))

    print("Finished Comparison 3\n")

    # Create a DataFrame to store results
    comparisons_df = pd.DataFrame(comparisons, columns=[
        "ComparisonID", "Study1", "Study2", "AvgOddsRatio", "ComparisonType"
    ])

    print("All comparisons completed. Returning results as a DataFrame.\n")

    return comparisons_df






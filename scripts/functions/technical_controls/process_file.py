#!/home/nelazzabi/anaconda3/envs/homl3/bin/python
import os
import scanpy as sc
import anndata as ad
import argparse
import sys

# Import functions from their respective modules
sys.path.append("/home/nelazzabi/rewiring/scripts")
from functions.technical_controls.shuffle_patients import filter_and_shuffle_patients
from functions.technical_controls.compute_cor import compute_gene_correlation
from functions.technical_controls.saveCorDict import cor_mtx_to_edgeList

def process_file(file_path, categories, n_shuffles, output_dir, correlation_type):
    """
    Process a single file to compute gene correlations, shuffle patients, and save the results.
    
    Arguments:
    - file_path (str): The path to the input file.
    - categories (set): The categories to filter by (e.g., AD, CTL).
    - n_shuffles (int): The number of shuffles to perform.
    - output_dir (str): Directory to save the results.
    - correlation_type (str): Type of correlation to compute (e.g., "pearson", "spearman").
    
    Returns:
    - None: The function processes and saves the result for a single file.
    """
    print(f"\nProcessing file: {file_path}")

    try:
        adata = sc.read_h5ad(file_path)
        available_categories = set(adata.obs["disease"].unique())

        if not categories.issubset(available_categories):
            print(f"Skipping {file_path}: Required categories {categories} not found ({available_categories})")
            return

        aggregated_matrix = {category: None for category in categories}
        gene_names = adata.var_names.tolist()

        # Shuffle patients and compute gene correlations for each shuffle
        for shuffle_idx in range(n_shuffles):
            print(f"  Shuffle {shuffle_idx + 1}/{n_shuffles}")
            shuffle_N = filter_and_shuffle_patients(adata, categories)

            for category in categories:
                if category in available_categories:
                    print(f"    Computing correlation for category: {category}")
                    shuffle_N_coExpr_matrix = compute_gene_correlation(
                        adata=shuffle_N,
                        correlation_type=correlation_type,
                        replace_nans=True,
                        min_cells_threshold=20,
                        category=category
                    )

                    if aggregated_matrix[category] is None:
                        aggregated_matrix[category] = shuffle_N_coExpr_matrix 
                    else:
                        aggregated_matrix[category] += shuffle_N_coExpr_matrix

        # Normalize the aggregated matrix by the number of shuffles
        for category in categories:
            aggregated_matrix[category] /= n_shuffles  # Normalize by number of shuffles

        # Save the results
        print(f"Saving aggregated matrix for {file_path}")
        edge_list_dir2 = os.path.join(output_dir, "edge-list", "all-matrix")
        os.makedirs(edge_list_dir2, exist_ok=True)

        for category in categories:
            agg_adata = ad.AnnData(X=aggregated_matrix[category])
            agg_adata.var_names = gene_names

            cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2)

        print(f"Successfully processed and saved results for {file_path}")

    except Exception as e:
        print(f"Skipping {file_path} due to error: {e}")


# Add an entry point to allow this script to be run as a standalone program
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a single file with specific parameters.")
    parser.add_argument('--file_path', required=True, help="Path to the input file")
    parser.add_argument('--categories', required=True, help="Categories to filter by, e.g., 'AD,CTL'")
    parser.add_argument('--n_shuffles', required=True, type=int, help="Number of shuffles")
    parser.add_argument('--output_dir', required=True, help="Directory to save results")
    parser.add_argument('--correlation_type', required=True, help="Type of correlation to use")

    args = parser.parse_args()

    # Convert the categories string to a set (if needed)
    categories = set(args.categories.split(','))
    process_file(args.file_path, categories, args.n_shuffles, args.output_dir, args.correlation_type)
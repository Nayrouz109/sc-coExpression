import os
import scanpy as sc
import anndata as ad
import glob
import sys

sys.path.append("/home/nelazzabi/rewiring/scripts")

# Import functions from their respective modules
from functions.technical_controls.shuffle_patients import filter_and_shuffle_patients
from functions.technical_controls.compute_cor import compute_gene_correlation
from functions.technical_controls.saveCorDict import cor_mtx_to_edgeList



def process_files(input_dir, cell_type, categories, n_shuffles, output_dir, correlation_type):
    print(f"Processing files in {input_dir} for cell type: {cell_type}")

    for file_path in glob.glob(os.path.join(input_dir, f"{cell_type}*.h5ad")):
        print(f"\nLoading file: {file_path}")

        try:
            adata = sc.read_h5ad(file_path)
            available_categories = set(adata.obs["disease"].unique())

            if not categories.issubset(available_categories):
                print(f"Skipping {file_path}: Required categories {categories} not found ({available_categories})")
                continue

            aggregated_matrix = {category: None for category in categories}
            gene_names = adata.var_names.tolist()

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

            for category in categories:
                aggregated_matrix[category] /= n_shuffles  # Normalize by number of shuffles

            print(f"Saving aggregated matrix for {file_path}")

            # Save aggregated matrix
            edge_list_dir2 = os.path.join(output_dir, "edge-list", "all-matrix")
            os.makedirs(edge_list_dir2, exist_ok=True)

            for category in categories:
                agg_adata = ad.AnnData(X=aggregated_matrix[category])
                agg_adata.var_names = gene_names
                
                cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2)

            print(f"Successfully processed and saved results for {file_path}")
        
        except Exception as e:
            print(f"Skipping {file_path} due to error: {e}")
            continue



#test = process_files(
#    input_dir= "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/00testing",
#    cell_type="Mic",  # Assuming “Mic” is the correct cell type
#    categories={"AD", "CTL"},  # Categories should be a set
#    n_shuffles=2,
#    output_dir= "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/00testing/out" ,
#    correlation_type= "pearson"  # Assuming this should be a string
#)
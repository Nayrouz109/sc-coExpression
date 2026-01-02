
import h5py
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from pathlib import Path
import json
import anndata as ad
import pickle
import os 

def cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2 ):

    # 1.  data into an edgeList
    df = pd.DataFrame(agg_adata.X, index=agg_adata.var_names, columns=agg_adata.var_names)
    edge_list = df.stack().reset_index()
    #edge_list.columns = ["source", "target", "weight"]
    edge_list.columns = ["geneA", "geneB", "rank_norm_corr"]
    
    # Remove self-loops & multiIndex
    edge_list_clean = edge_list[edge_list["geneA"] != edge_list["geneB"]]
    edge_list_clean.set_index(['geneA', 'geneB'], inplace=True)
    
    # Count the number of targets per source
    target_counts = edge_list_clean.groupby(level=0).size()
    target_counts_df = target_counts.reset_index()
    target_counts_df.columns = ["geneA", "num_targets"]
    
    # Sanity check
    if (target_counts_df["num_targets"] == (agg_adata.shape[0]-1)).all():
    #if (target_counts_df["num_targets"] == 10907).all():
        print(f"Sanity check passed")
        
        # Save edge list 
        output_file = edge_list_dir2 / f"{file_path.stem}_edgeList.csv"
        edge_list_clean.to_csv(output_file)
        print(f"Saved edge list to: {output_file}\n")
    else:
        print(f"Sanity check failed. Skipping file.\n")


def correlation_matrix_to_edgelist(corr_matrix: pd.DataFrame):

    upper_triangle_mask = np.triu(np.ones(corr_matrix.shape), k=1) ## mask the upper triangle excluding the diagonal
    
    # Apply the mask to the correlation matrix to get the non-redundant values
    upper_triangle_corr = corr_matrix.where(upper_triangle_mask == 1)
    
    # Unstack the upper triangle into a series, which will automatically drop NaN values
    edge_list = upper_triangle_corr.stack().reset_index()
    
    # Rename columns to match the desired edge list format
    edge_list.columns = ['geneA', 'geneB', 'rank_normalized_correlation']
    
    # Optionally, you can set geneA and geneB as the index
    edge_list.set_index(['geneA', 'geneB'], inplace=True)
    
    return edge_list

def save_to_hdf5(cor, output_dir, file_path):
    """
    Save correlation data, including edge lists and aggregated sparse matrices, to specified directories.
    
    Args:
        cor (dict): The correlation dictionary containing matrices, gene names, and data summary.
        output_dir (Path): The base output directory to save the data.
        file_path (Path): The input file path (used for naming output files).
    """
    # Create required directories
    agg_mtx_dir = output_dir / "agg-mtx"
    edge_list_dir1 = output_dir / "edge-list" / "upTriang"
    edge_list_dir2 = output_dir / "edge-list" / "all-matrix"
    summary_dir = output_dir / "cor-summary"

    agg_mtx_dir.mkdir(parents=True, exist_ok=True)
    edge_list_dir1.mkdir(parents=True, exist_ok=True)
    edge_list_dir2.mkdir(parents=True, exist_ok=True)
    summary_dir.mkdir(parents=True, exist_ok=True)

    # File names
    agg_file = agg_mtx_dir / f"{file_path.stem}_corAggSparse.h5ad"
    edge_list_file = edge_list_dir1 / f"{file_path.stem}_corEdgeLists.csv"
    #edge_list_file = edge_list_dir / f"{file_path.stem}_corEdgeLists.csv.gz"
    df_summary_file = summary_dir / f"{file_path.stem}_corSummary.csv"

    #1. Save AnnData object
    aggregated_matrix = cor["aggregated_rank_normalized_matrix"]

    gene_names = cor["gene_names"]
    agg_adata = ad.AnnData(X=aggregated_matrix)
    agg_adata.var_names = gene_names
    agg_adata.write_h5ad(agg_file, compression="gzip")
    print(f"Saved aggregated matrix and gene names to {agg_file}")

    #2. Save edgeList
    aggregated_matrixEdge = cor["edgeList"] 
    edge_list= correlation_matrix_to_edgelist(aggregated_matrixEdge)
    edge_list.to_csv(edge_list_file)
    print("saved edge list")

    #2.2 Save edgeList
    cor_mtx_to_edgeList(file_path, agg_adata, edge_list_dir2)

    #3. Save data summary 
    #cor["data_summary"].to_csv(df_summary_file, index=False)
    #print(f"Saved data summary to {df_summary_file}")
    
















    #edgeList = cor["corEdgeList_per_patientID"] 
    #edgeList.to_csv(edge_list_file, index=False, compression='gzip')

    # Save raw edge list as JSON
    #with open(correlation_edge_list_file, 'wb') as f:
        #pickle.dump(cor["corEdgeList_per_patientID"], f)
    #print(f"Saved raw edge list to {correlation_edge_list_file}")
    
    #with open(rank_edge_list_file, 'wb') as f:
        #pickle.dump(cor["rnkEdgeList_per_patientID"], f)
    #print(f"Saved ranked edge list to {rank_edge_list_file}")

    #with open(correlation_edge_list_file, 'w') as f:
        #json.dump(convert_tuple_keys(cor["corEdgeList_per_patientID"]), f, indent=4)
    # Save ranked edge list as JSON
    #with open(rank_edge_list_file, 'w') as f:
        #json.dump(convert_tuple_keys(cor["rnkEdgeList_per_patientID"]), f, indent=4)
    


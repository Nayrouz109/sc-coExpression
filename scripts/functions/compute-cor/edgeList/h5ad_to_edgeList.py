
import os
import anndata as ad
import pandas as pd

def process_h5ad_file(input_path, output_path):
    print(f"Processing: {input_path}")
    

    adata = ad.read_h5ad(input_path)
    
    # 1.  data into an edgeList
    df = pd.DataFrame(adata.X, index=adata.var_names, columns=adata.var_names)
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
    if (target_counts_df["num_targets"] == 10907).all():
        print(f"Sanity check passed for: {input_path}")
        
        # Save edge list
        output_file = os.path.join(output_path, os.path.basename(input_path).replace("_corAggSparse.h5ad", "_edgeList.csv"))
        edge_list_clean.to_csv(output_file)
        print(f"Saved edge list to: {output_file}\n")
    else:
        print(f"Sanity check failed for: {input_path}. Skipping file.\n")
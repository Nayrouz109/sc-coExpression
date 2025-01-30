


df5 = cor5["AD"]["aggregated_rank_normalized_matrix"]

import pandas as pd
import numpy as np

def correlation_matrix_to_edgelist(corr_matrix: pd.DataFrame):

    upper_triangle_mask = np.triu(np.ones(corr_matrix.shape), k=1) ## mask the upper triangle excluding the diagonal
    
    # Apply the mask to the correlation matrix to get the non-redundant values
    upper_triangle_corr = corr_matrix.where(upper_triangle_mask == 1)
    
    # Unstack the upper triangle into a series, which will automatically drop NaN values
    edge_list = upper_triangle_corr.stack().reset_index()
    
    # Rename columns to match the desired edge list format
    edge_list.columns = ['geneA', 'geneB', 'rank_normalized_correlation']
    
    # Optionally, you can set geneA as the index
    #edge_list.set_index('geneA', inplace=True)
    #edge_list = edge_list.set_index(['geneA', 'geneB']).reset_index()
    edge_list.set_index(['geneA', 'geneB'], inplace=True)
    
    return edge_list

# Example usage:
# Assuming `corr_matrix` is your DataFrame of shape (10908, 10908)
edge_list = correlation_matrix_to_edgelist(df5)

### to access a gene 
edge_list.loc["FAM161A"]

### save the edge list 
edgeList.to_csv("/home/nelazzabi/edgeList.csv") 
loaded_edge_list = pd.read_csv('/home/nelazzabi/edgeList.csv', index_col=[0, 1])


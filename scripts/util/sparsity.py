

import scipy.sparse

def compute_sparsity(sparse_matrix):
    """
    Computes the percentage of sparsity in a given sparse matrix.
    
    Parameters:
        sparse_matrix (scipy.sparse.spmatrix): A sparse matrix.
        
    Returns:
        float: The percentage of sparsity in the matrix.
    """
    if not scipy.sparse.issparse(sparse_matrix):
        raise ValueError("Input must be a sparse matrix.")
    
    total_elements = sparse_matrix.shape[0] * sparse_matrix.shape[1]
    non_zero_elements = sparse_matrix.nnz  # Number of non-zero elements
    
    sparsity = (1 - (non_zero_elements / total_elements)) * 100  # Percentage of zeros
    return sparsity


#### example run 
from scipy.sparse import csr_matrix

# Example sparse matrix
sparse_mat = csr_matrix([[0, 1, 0], [0, 0, 0], [2, 0, 0]])

# Compute sparsity percentage
sparsity_percentage = compute_sparsity(sparse_mat)
print(f"Sparsity: {sparsity_percentage:.2f}%")







import numpy as np
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt

### to get the data, refer back to testing.py script 
def sparse_corr(A):
    """C
    Jan 27th, 2025 
    """
    A = A.T  # Transpose for gene-wise computation
    mean = A.mean(axis=1)
    A = csr_matrix(A - mean)
    norm = np.sqrt(A.multiply(A).sum(axis=1))
    A = A.multiply(1 / norm)
    return A @ A.T  # Sparse matrix multiplication


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

def sparse_corrV2(A):
    """Compute Pearson correlation for a sparse matrix.
    
    Arguments:
        A {scipy.sparse matrix} -- Sparse matrix of gene expression values
    Returns:
        np.matrix -- Correlation matrix
    """
    N = A.shape[0] # Calculate the number of rows (cells)

    C = ((A.T *A -(sum(A).T*sum(A)/N))/(N-1)).todense() # Calculate the covariance matrix
    V = np.sqrt(np.mat(np.diag(C)).T * np.mat(np.diag(C))) # Compute the standard deviation matrix

    COR = np.divide(C, V + 1e-119) # Compute the correlation matrix

    # Find genes (columns in A) with a sum of 0
    zero_sum_columns = (np.sum(A, axis=0) == 0)  # Sum across rows (cells)

    # Convert the boolean array to an integer index (this will ensure it can be used for indexing)
    zero_sum_columns = np.where(zero_sum_columns)[0]  # Get the indices of zero-sum columns

    # Set correlations involving these genes (rows and columns in COR) to NaN
    COR[zero_sum_columns, :] = np.nan  # Set rows corresponding to zero-sum genes to NaN
    COR[:, zero_sum_columns] = np.nan  # Set columns corresponding to zero-sum genes to NaN

    return COR

dense_matrix = subset_data.toarray()


corr_matrix1 = sparse_corr(subset_data)
corr_matrix2 = sparse_corrV2(subset_data)
corr_matrix3 = np.corrcoef(dense_matrix, rowvar=False)


### check if the cor matrices has NAs ===========================
np.fill_diagonal(corr_matrix1, 1)
np.fill_diagonal(corr_matrix2, 1)
np.fill_diagonal(corr_matrix3, 1)

np.sum(np.isnan(corr_matrix1))
np.sum(np.isnan(corr_matrix2))
np.sum(np.isnan(corr_matrix3))

### counts 0s 
gene_sums = subset_data.sum(axis=0).A1  # Sum of each gene across all cells
zero_sum_genes = (gene_sums == 0).sum()

### plot his distributions 
plt.hist(np.ravel(corr_matrix1), bins=30, alpha=0.5, label='Matrix 1')
plt.hist(np.ravel(corr_matrix2), bins=30, alpha=0.5, label='Matrix 3')

### rank the matrix ====================================
np.fill_diagonal(corr_matrix3, 1)


corr_matrix3R =  rank(corr_matrix3)
np.sum(np.isnan(corr_matrix3R))
plt.hist(np.ravel(corr_matrix3R), bins=30)
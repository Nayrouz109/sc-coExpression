



import numpy as np    
from scipy import sparse

def sparse_corr(A):
    N = A.shape[0]
    C=((A.T*A -(sum(A).T*sum(A)/N))/(N-1)).todense()
    V=np.sqrt(np.mat(np.diag(C)).T*np.mat(np.diag(C)))
    COR = np.divide(C,V+1e-119)
    return COR



A = sparse.rand(1000000, 100, density=0.1, format='csr')
test1 = sparse_corr(A)
test2 = np.corrcoef(A.T.toarray())

## source: https://stackoverflow.com/questions/19231268/correlation-coefficients-for-sparse-matrix-in-python



M_done = compute_gene_correlation(M, 
                                  correlation_type='pearson', 
                                  replace_nans=True, 
                                  min_cells_threshold=20, 
                                  max_workers=7)

M1_done = process_patient(patient_id, 
                          M, 
                          min_cells_threshold=20, 
                          correlation_type='pearson', 
                          replace_nans=True)





### working with annData ****************************
# one plot 
plt.hist(CT.X.data, bins=50)  # You can adjust the number of bins
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram of Non-Zero Values in Sparse Matrix')
plt.show()
# two plots 
plt.figure(figsize=(8, 6))
plt.hist(AD.X.data, bins=30, alpha=0.5, label='Matrix 1', color='blue')
plt.hist(CT.X.data, bins=30, alpha=0.5, label='Matrix 2', color='orange')

### working with scipy.sparse.csr_matrix ****************************
plt.hist(corr_matrix.data)

### assuming that you have a dictionary and you want to plot a histogram 
### distribution of its correlation values 

import matplotlib.pyplot as plt
import numpy as np

correlation_matrix_AD = cor["AD"]["aggregated_rank_normalized_matrix"] 
correlation_matrix_CT= cor["CTL"]["aggregated_rank_normalized_matrix"] 

# Convert the sparse matrices to dense format
dense_correlation_matrix_AD = correlation_matrix_AD.toarray()
dense_correlation_matrix_CT = correlation_matrix_CT.toarray()

# Flatten the dense matrices to get all correlation values
correlation_values_AD = dense_correlation_matrix_AD.flatten()
correlation_values_CT = dense_correlation_matrix_CT.flatten()

# Plot the histograms on top of each other
plt.figure(figsize=(8, 6))
plt.hist(correlation_values_AD, bins=50, alpha=0.5, label='AD', color='blue', edgecolor='black')
plt.hist(correlation_values_CT, bins=50, alpha=0.5, label='CTL', color='green', edgecolor='black')
plt.title('Histogram of Correlation Values: AD vs CTL')
plt.xlabel('Correlation Value')
plt.ylabel('Frequency')
plt.legend()
plt.grid(True)
plt.show()


#### working with panda **************************************
all_values = df.values.flatten() 

#### **************************************


import numpy as np
import matplotlib.pyplot as plt




### numpy.matrix *************************
corr_flat = np.ravel(norm)

# Plot the histogram of the correlation matrix
plt.figure(figsize=(8, 6))
plt.hist(corr_flat, bins=50, color='blue', alpha=0.7)
plt.title('Histogram of Correlation Matrix')
plt.xlabel('Correlation Values')
plt.ylabel('Frequency')
plt.show()

# if you want to plot the dist of 2 matrices 
corr_flat1 = np.ravel(test1.X)
corr_flat2 = np.ravel(test2.X)
corr_flat3 = np.ravel(test3.X)

# Plot the histogram of the correlation matrix
plt.figure(figsize=(8, 6))
plt.hist(corr_flat1, bins=50, color='blue', alpha=0.5, label='Dataset 1')
plt.hist(corr_flat2, bins=50, color='orange', alpha=0.5, label='Dataset 2')
plt.title('Histograms of Correlation Matrices')
plt.xlabel('Correlation Values')
plt.ylabel('Frequency')
plt.legend()  # Add a legend to differentiate the datasets
plt.show()











# Assuming M1_done[1]["correlation_matrix"] and M1_done[1]["rank_normalized_matrix"] are matrices
corr_matrix = M1_done[1]["correlation_matrix"]
rank_norm_matrix = M1_done[1]["rank_normalized_matrix"]

# Flatten the matrices to 1D arrays
corr_flat = np.ravel(corr_matrix)
rank_norm_flat = np.ravel(rank_norm_matrix)

# Create a figure with two subplots for the histograms
plt.figure(figsize=(12, 6))

# Plot histogram for correlation matrix
plt.subplot(1, 2, 1)  # First subplot
plt.hist(corr_flat, bins=50, color='blue', alpha=0.7)
plt.title('Histogram of Correlation Matrix')
plt.xlabel('Correlation Values')
plt.ylabel('Frequency')

# Plot histogram for rank normalized matrix
plt.subplot(1, 2, 2)  # Second subplot
plt.hist(rank_norm_flat, bins=50, color='green', alpha=0.7)
plt.title('Histogram of Rank Normalized Matrix')
plt.xlabel('Rank Normalized Values')
plt.ylabel('Frequency')

# Show the plot
plt.tight_layout()
plt.show()




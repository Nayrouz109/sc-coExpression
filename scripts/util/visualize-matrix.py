

###1. matplotlib

import matplotlib.pyplot as plt
import numpy as np

# Example matrix
matrix = np.random.rand(10, 10)  # Replace with your matrix

plt.imshow(matrix, cmap='viridis', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('Matrix Heatmap')
plt.xlabel('Column Index')
plt.ylabel('Row Index')
plt.show()

###2. seaborn (w/ enhanced features)
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# Example matrix
matrix = np.random.rand(10, 10)  # Replace with your matrix

sns.heatmap(matrix, annot=True, fmt=".2f", cmap="coolwarm", cbar_kws={'label': 'Value'})
plt.title('Matrix Heatmap')
plt.xlabel('Column Index')
plt.ylabel('Row Index')
plt.show()
 

### 3. 3D w/ matplotlib 
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import matplotlib.pyplot as plt

# Example matrix
matrix = np.random.rand(10, 10)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(range(matrix.shape[1]), range(matrix.shape[0]))
ax.plot_surface(X, Y, matrix, cmap='viridis')
ax.set_title('Matrix in 3D')
plt.show()



### 4. interactive w/ plotly 
import plotly.express as px
import numpy as np

# Example matrix
matrix = np.random.rand(10, 10)

fig = px.imshow(matrix, color_continuous_scale='viridis', labels={'color': 'Value'})
fig.update_layout(title='Matrix Heatmap', xaxis_title='Column Index', yaxis_title='Row Index')
fig.show()



### you'd want to make sure matrix is sparse 
subset_data = subset_data.toarray()
### also look at fewer cells 

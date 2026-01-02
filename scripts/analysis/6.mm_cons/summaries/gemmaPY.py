#!/usr/bin/env python
# coding: utf-8

# In[5]:


import gemmapy
from getpass import getpass


# In[66]:

# add username to see private datasets
# just uncomment that part of the line and add username (no need to add password)
api = gemmapy.GemmaPy(path='https://staging-gemma.msl.ubc.ca/rest/v2')#, auth=(input('Username: '), getpass('Password: ')))


# In[67]:

# get all datasets with category='assayType' and value = 'single nuclei sequencing' or 'single cell sequencing'
# be careful, because some of them are not actually single cell datasets according to guillaume (this seems wrong but he doesn't care)
# when you try to get single cell dimensions, the command will fail for those datasets

data = api.get_all_pages(api.get_datasets,filter='allCharacteristics.categoryUri = http://purl.obolibrary.org/obo/OBI_0000070 and allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/OBI_0002631, http://purl.obolibrary.org/obo/OBI_0003109)', sort='-id')
# returns a data frame that you can filter by taxon

# In[68]:

# public experiment ocunts by taxom
data.taxon_name.value_counts()


# In[69]:

# get sample count histogram
data.experiment_sample_count.hist(bins=70)


# In[70]:


data.head(20)

filtered_df = data[data['taxon_scientific'] == 'Mus musculus']
num_unique_experiments = df[df['experiment_sample_count'] < 10]['experiment_short_name'].nunique()

# In[71]:

# example of a single cell dataset
# you could iterate through all mouse datasets and insert a try except when fetching single cell dimensions
api.raw.get_dataset_single_cell_dimension('GSE296027').data.number_of_cells

# Define a function that gets the number of cells for a dataset
def get_cell_count(experiment_name):
    return api.raw.get_dataset_single_cell_dimension(experiment_name).data.number_of_cells

# Apply the function to the 'experiment_short_name' column and create a new column
filtered_df['total_cell_count'] = filtered_df['experiment_short_name'].apply(get_cell_count)

#GSE269482


# In[72]:

# for cell type assignment in a given dataset, get number of cell types
for cta in api.raw.get_dataset_single_cell_dimension('GSE296027').data.cell_type_assignments:
    print(cta.cell_types)


# In[74]:

# get sample factor values for a given dataset
api.get_dataset_samples('GSE296027').sample_factor_values[0]


# In[75]:

# get sample characteristics for a given dataset
api.get_dataset_samples('GSE296027').sample_characteristics[0]


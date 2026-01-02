

import gemmapy
from getpass import getpass


api = gemmapy.GemmaPy(path='https://staging-gemma.msl.ubc.ca/rest/v2')#, auth=(input('Username: '), getpass('Password: ')))


data = api.get_all_pages(api.get_datasets,filter='allCharacteristics.categoryUri = http://purl.obolibrary.org/obo/OBI_0000070 and allCharacteristics.valueUri in (http://purl.obolibrary.org/obo/OBI_0002631, http://purl.obolibrary.org/obo/OBI_0003109)', sort='-id')
data.taxon_name.value_counts()

### sample count - other unique values count 
mouse_data = data[data["taxon_name"].str.contains("mouse", na=False)]
total_samples = mouse_data["experiment_sample_count"].sum() # total samples 
mouse_data.experiment_sample_count.hist(bins=70)            # total samples plot 
mouse_data[mouse_data['experiment_sample_count'] >20]['experiment_short_name'].nunique()

unique_counts = {col: mouse_data[col].nunique() for col in ["experiment_short_name", "experiment_ID", "experiment_accession"]}


### single cells count - add column 
from gemmapy.sdk.rest import ApiException

def get_cell_count(experiment_name):
    try:
        return api.raw.get_dataset_single_cell_dimension(experiment_name).data.number_of_cells
    except ApiException:
        return None  # or 0 if you'd rather treat missing as zero

mouse_data["total_cell_count"] = mouse_data["experiment_short_name"].apply(get_cell_count)

#def get_cell_count(experiment_name):
#    return api.raw.get_dataset_single_cell_dimension(experiment_name).data.number_of_cells
#mouse_data['total_cell_count'] = mouse_data['experiment_short_name'].apply(get_cell_count)

mouse_data.total_cell_count.hist(bins=70) # plot cells count 
mouse_data[mouse_data['total_cell_count'] >20]['experiment_short_name'].nunique() 

### cell type annotations 
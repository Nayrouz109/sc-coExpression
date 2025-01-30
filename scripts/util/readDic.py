

### using json module 
import json

gene_dict_path = '/home/nelazzabi/rewiring/data/genes/filtered_genes.json'
with open(gene_dict_path, "r") as fp: 
    gene_dict = json.load(fp) 


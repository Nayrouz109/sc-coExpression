
import h5py
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from pathlib import Path
import json
import anndata as ad
import pickle



def save_to_csv(cor, output_dir, file_path):

    edge_list_dir = output_dir / "edge-list" 
    edge_list_dir.mkdir(parents=True, exist_ok=True)

    #edge_list_file = edge_list_dir / f"{file_path.stem}_corEdgeLists.csv.gz"
    edge_list_file = edge_list_dir / f"{file_path.stem}_corEdgeLists.csv"

    edgeList = cor["corEdgeList_per_patientID"] 
    #edgeList.to_csv(edge_list_file, index=False, compression='gzip')
    #edgeList.to_csv(edge_list_file, index=False)
    edgeList.to_csv(edge_list_file) 




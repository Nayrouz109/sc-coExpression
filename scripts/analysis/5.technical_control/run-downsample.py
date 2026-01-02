#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

import os
import sys


# Define parameters
input_dir = "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease"
output_dir_samples = os.path.join(input_dir, "02downSampled/dsp-n")
output_dir_cells = os.path.join(input_dir, "02downSampled/AD-Ctl-smpls-clls")
cell_type = "Mic"
conditions = ["AD", "CTL"]
downsample_cells=False
iterations=50 

sys.path.append("/home/nelazzabi/rewiring/scripts")
#from functions.technical_controls.downSample import process_files
#process_files(input_dir, output_dir_samples, output_dir_cells, cell_type, conditions)


from functions.technical_controls.NdownSample import process_files

# Run the process
process_files(input_dir, 
              output_dir_samples, 
              output_dir_cells, 
              cell_type, 
              conditions, 
              downsample_cells,
              iterations)

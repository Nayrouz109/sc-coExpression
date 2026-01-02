#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

import os
import glob
import sys 
import subprocess
from pathlib import Path

sys.path.append("/home/nelazzabi/rewiring/scripts")
from functions.technical_controls.process_file import process_file

def submit_slurm_job(file_path, categories, n_shuffles, output_dir, correlation_type):
    """
    Generate and submit a SLURM job to process a single file using the `process_file` function.
    """
    #slurm_script = f"""
    job_name = f"cor_{Path(file_path).stem}"
    categories_str = ','.join(categories) 
    slurm_script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={output_dir}/logs/{job_name}.out
#SBATCH --error={output_dir}/logs/{job_name}.err
#SBATCH --cpus-per-task=20


# Run the Python script for this file
python /home/nelazzabi/rewiring/scripts/functions/technical_controls/process_file.py --file_path {file_path}  --categories {categories_str} --n_shuffles {n_shuffles} --output_dir {output_dir} --correlation_type {correlation_type}
""" 
    
    # Write SLURM script to a file
    slurm_filename = f"process_file_{os.path.basename(file_path)}.slurm"
    with open(slurm_filename, 'w') as slurm_file:
        slurm_file.write(slurm_script)

    # Submit the SLURM job using sbatch
    subprocess.run(['sbatch', slurm_filename])

    # Optionally, remove the SLURM script after submission
    #os.remove(slurm_filename)

def process_files(input_dir, cell_type, categories, n_shuffles, output_dir, correlation_type):
    """
    Loop through all files in input_dir and submit SLURM jobs for each file.

    Arguments:
    - input_dir: Directory containing input files.
    - cell_type: Cell type to process.
    - categories: Categories to filter by.
    - n_shuffles: Number of shuffles to perform.
    - output_dir: Directory to save results.
    - correlation_type: Correlation type (e.g., Pearson).
    """
    print(f"Processing files in {input_dir} for cell type: {cell_type}")
    
    # Loop over all .h5ad files and submit SLURM jobs
    for file_path in glob.glob(os.path.join(input_dir, f"{cell_type}*.h5ad")):
        submit_slurm_job(file_path, categories, n_shuffles, output_dir, correlation_type)

if __name__ == "__main__":
    #input_dir = "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease"
    input_dir = "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/00testing" 
    cell_type = "Mic"  
    categories = 'AD','CTL'
    n_shuffles = 2  
    output_dir = "/cosmos/data/project-data/NW-rewiring/coExpr/n-shuffled"
    correlation_type = "pearson" 
    
    # Call the process_files function to start the SLURM submission
    process_files(input_dir, cell_type, categories, n_shuffles, output_dir, correlation_type)

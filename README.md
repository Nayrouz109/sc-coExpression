# gene co-expression "network" Rewiring

#### project purpose and goals 
1- NW in Disease (AD and SZC)
2- NW in Cell types ("normal" samples)

The project does the following: 
1- pipeline to compute gene co-expression matrices 
2- pipeline to compute differential gene co-expression (absolute difference)
3- pipeline to compute metrics for differential co-expression 
3.1 co-expression conservation score 
3.2 FDR 
4- pipeline to compute differential gene expression 

### where to find data and scripts 
# data:   /cosmos/data/project-data/NW-rewiring/data
# result: /cosmos/data/project-data/NW-rewiring/coExpr
# script: /home/nelazzabi/rewiring/scripts


# Installation 
1. clone the repo 
git clone https://github.com/yourusername/bioinformatics-project.git
2. Set up the Conda environment:
conda env create -f envPython.yml conda activate bioinformatics-env
3. run the main analysis cor script - ex AD 
python run_compute_cor-disease.py \
    --input_dir  /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease \
    --output_dir /cosmos/data/project-data/NW-rewiring/coExpr/disease \
    --correlation_type pearson \
    --replace_nans True \
    --analysis_type disease \
    --categories AD CTL
4. 


# Project structure 
- data\
- scripts\


# Pipeline notes 
the output from the parallel computation 
            subset_results = Parallel(n_jobs=max_workers)(
                delayed(process_patient)(patient_id, subset, min_cells_threshold, correlation_type, replace_nans)
                for patient_id in patient_ids
            )

subset_results, is expected to be a collection (like a list or a generator) of tuples. Each tuple contains two elements:
	patient_id: This is the identifier for a patient.
    	result: This is a dictionary or object containing the processed data for that patient

subset_results is a list 
subset_results[0]    # Access the first tuple - data for the first patient 
subset_results[0][0] # Access the identifier ('82317494')
subset_results[0][1] # Access the dictionary for the patient
subset_results[0][1]["rank_normalized_matrix"] 

# Contact 
for questions, contact nayrouz@student.ubc.ca


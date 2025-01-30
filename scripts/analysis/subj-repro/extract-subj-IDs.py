
import pandas as pd
import os
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path


### this is to extract patient IDs to work with 

#   1. only work with control samples to avoide inter-individual disease variability 
#   2. i chose to work with PFC samples that are sequenced in 2 studies 
#   3. but also sample samples are seq in 6 brain regions in one of the studies 


p_meta = pd.read_csv("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData.csv")
new_p = p_meta[(p_meta["study"] != "ROSMAP-microGlia") & (p_meta["diagnosis"] == "no")]
data = new_p

# Summarize occurrences of patientID
summary = (
    data.groupby('patientID')
    .agg(
        count=('patientID', 'size'),
        unique_brain_regions=('brainRegion', 'nunique'),
        unique_studies=('study', 'nunique')
    )
    .reset_index()
)

# Add a column to describe the context
summary['context'] = summary.apply(
    lambda row: (
        "Appears in multiple brain regions within the same study"
        if row['unique_studies'] == 1 and row['unique_brain_regions'] > 1
        else "Appears across multiple studies"
        if row['unique_studies'] > 1
        else "Appears in a single brain region within a single study"
    ),
    axis=1
)

summaryV1 = summary[summary["count"] > 1]
summaryV2 = summaryV1[summaryV1["count"] == 7]

patientIDs = summaryV2["patientID"].unique()

df = pd.DataFrame(patientIDs) 
df.columns = ["patientID"]
df.to_csv("/home/nelazzabi/rewiring/scripts/analysis/subj-repro/patientIDs.csv", header= True, index= False)

dfV2 = pd.read_csv("/home/nelazzabi/rewiring/scripts/analysis/subj-repro/patientIDs.csv")


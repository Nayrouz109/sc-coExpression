



import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Load the dataset
file_path = "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/all-patient-metaData.csv"
p_meta = pd.read_csv(file_path)

study_summary = (
    p_meta.groupby('study')
    .size()
    .reset_index(name='Number of Samples')
    .sort_values(by='Number of Samples', ascending=False)
)

### single cells ==============================================
sc_meta = pd.read_csv("/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag/summary.csv")
sc_meta_summary = (sc_meta.groupby('study name')['total single cells'].sum()
                   .reset_index().
                   sort_values(by= "total single cells", ascending=False))



# Load the dataset =================================================

summary = (
    p_meta.groupby(['study', 'diagnosis'])
    .size()
    .unstack(fill_value=0)
    .rename(columns={'yes': 'disease', 'no': 'control'})
)

# Add a total count column
summary['Total Samples'] = summary['disease'] + summary['control']

# Add the disease column (assuming all patients in the same study have the same disease)
summary['Disease'] = (
    p_meta.groupby('study')['disease']
    .first()  # Assumes all rows for the same study have the same disease value
)

# Order by total sample count (descending)
summary = summary.sort_values(by='Total Samples', ascending=False).reset_index()
summary.loc[summary['study'].isin(['MSBB-2024', 'HBCC-2024', 'RADC-2024']), 'Disease'] = 'mix'

summary = summary[summary["study"] != "ROSMAP-microGlia"]

total_summary = pd.merge(summary, 
                         sc_meta_summary, 
                         left_on= "study", 
                         right_on="study name", 
                         how= "inner") 

example_data = testF[testF["study name"].isin(["ROSMAP-DeJager-2024", "ROSMAP-Mathys-2023", "Kamath-NPH-2023", "ROSMAP-erosion-2023", "ROSMAP-Mathys-2019"])]
# Save the summary table
#summary_file = "/space/scratch/nairuz-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent/study_summary_disease_control.csv"
#summary.to_csv(summary_file, index=False)





### distribution of sex per diagnosis per study ==================================
p_meta['diagnosis'] = p_meta['diagnosis'].str.lower()
grouped_data = (
    p_meta.groupby(['study', 'diagnosis', 'sex'])
    .size()
    .reset_index(name='count')
)

# Pivot the data for easier plotting
pivot_data = grouped_data.pivot_table(
    index='study',
    columns=['diagnosis', 'sex'],
    values='count',
    fill_value=0
)

# Define variables for plotting
studies = pivot_data.index
diagnoses = ['yes', 'no']
sexes = ['male', 'female']
bar_width = 0.4
x = np.arange(len(studies))  # Positions for the study groups

# Define color-blind-friendly colors
male_color = "#2166AC"  # Dark blue
female_color = "#B2182B"  # Dark red

# Initialize the plot
fig, ax = plt.subplots(figsize=(12, 8))
ax.set_facecolor('white') # Remove grid and add a clean white background
fig.patch.set_facecolor('white')
ax.grid(False)

# Plot bars for each diagnosis
for i, diagnosis in enumerate(diagnoses):
    male_counts = pivot_data[(diagnosis, 'male')].values
    female_counts = pivot_data[(diagnosis, 'female')].values

    # Bottom of the stack (start with males)
    ax.bar(
        x + i * bar_width, 
        male_counts, 
        bar_width, 
        label=f'{diagnosis.capitalize()} - Male', 
        color=male_color
    )
    
    # Stack females on top of males
    ax.bar(
        x + i * bar_width, 
        female_counts, 
        bar_width, 
        bottom=male_counts, 
        label=f'{diagnosis.capitalize()} - Female', 
        color=female_color
    )

# Customize the plot
ax.set_xticks(x + bar_width / 2)
ax.set_xticklabels(studies, rotation=45, ha='right')
ax.set_xlabel("Study", fontsize=12)
ax.set_ylabel("Number of Patients", fontsize=12)
ax.set_title("Male and Female Counts by Diagnosis per Study (Stacked Bar)", fontsize=14)
ax.legend(title="Diagnosis and Sex", fontsize=10)

# Adjust spacing
plt.tight_layout()

# Save and show the plot
plt.savefig("stacked_bar_plot_grouped_by_study_clean.png", dpi=300)
plt.show()






# Set seaborn style
sns.set(style="whitegrid")

# Plot 1: Age distribution
plt.figure(figsize=(8, 6))
sns.histplot(p_meta['age'], bins=20, kde=True, color='skyblue')
plt.title("Age Distribution of Patients")
plt.xlabel("Age")
plt.ylabel("Frequency")
plt.savefig("age_distribution.png")
plt.show()

# Plot 2: Sex distribution
plt.figure(figsize=(6, 6))
sns.countplot(data=p_meta, x='sex', palette='pastel')
plt.title("Sex Distribution of Patients")
plt.xlabel("Sex")
plt.ylabel("Count")
plt.savefig("sex_distribution.png")
plt.show()

# Plot 3: Diagnosis and Disease distribution
plt.figure(figsize=(10, 6))
sns.countplot(data=p_meta, x='diagnosis', hue='disease', palette='muted')
plt.title("Diagnosis and Disease Distribution")
plt.xlabel("Diagnosis")
plt.ylabel("Count")
plt.legend(title="Disease")
plt.savefig("diagnosis_disease_distribution.png")
plt.show()

# Plot 4: Patients by Brain Region and Study
plt.figure(figsize=(12, 8))
sns.countplot(data=p_meta, x='brainRegion', hue='study', palette='deep')
plt.title("Patients by Brain Region and Study")
plt.xlabel("Brain Region")
plt.ylabel("Count")
plt.legend(title="Study")
plt.xticks(rotation=45, ha='right')
plt.savefig("brain_region_study_distribution.png")
plt.show()


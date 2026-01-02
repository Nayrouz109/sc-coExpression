#!/home/nelazzabi/anaconda3/envs/homl3/bin/python

import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd

### mainting true ratio of AD and CTL per study 

def parse_args():
    parser = argparse.ArgumentParser(description="Shuffle disease annotation for single cells from each patient.")
    parser.add_argument("--input_dir", type=str, required=True, help="Directory containing .h5ad files")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save shuffled .h5ad files")
    parser.add_argument("--categories", type=str, required=True, help="Comma-separated disease categories to filter by (e.g., AD,CTL)")
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    categories = set(args.categories.split(","))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith(".h5ad"):
            file_path = os.path.join(input_dir, filename)
            print(f"\nProcessing file: {filename}\n")

            # Load .h5ad file
            adata = sc.read_h5ad(file_path)
            
            if "patientID" not in adata.obs.columns or "disease" not in adata.obs.columns:
                print(f"Skipping {filename}: Missing 'patientID' or 'disease' column in .obs")
                continue

            unique_diseases = set(adata.obs["disease"].unique())
            # Ensure both categories exist
            if not categories.issubset(unique_diseases):
                print(f"Skipping {filename}: Required categories {categories} not found in data ({unique_diseases})")
                continue

            # Filter patient IDs based on disease categories
            valid_patients = adata.obs[adata.obs["disease"].isin(categories)]["patientID"].unique()
            filtered_data = adata[adata.obs["patientID"].isin(valid_patients)]

            # Shuffle patients into two random groups
            np.random.seed(42)  # Ensure reproducibility
            all_patients = list(valid_patients)
            np.random.shuffle(all_patients)
            # Calculate patient numbers for each disease category
            num_ctl = np.sum(adata.obs.drop_duplicates(subset='patientID')["disease"] == "CTL")
            num_ad = np.sum(adata.obs.drop_duplicates(subset='patientID')["disease"] == "AD")
            grp1_patients = set(all_patients[:num_ctl])
            grp2_patients = set(all_patients[num_ctl:num_ctl+num_ad])

            #num_ctl = np.sum(filtered_data.obs["disease"] == "CTL")
            #num_ad = np.sum(filtered_data.obs["disease"] == "AD")

            # Group patients maintaining the true ratio
            #grp1_patients = set(filtered_data.obs["patientID"][:num_ctl])
            #grp2_patients = set(filtered_data.obs["patientID"][num_ctl:num_ctl+num_ad])

            # Assign random group labels
            filtered_data.obs["random_group"] = filtered_data.obs["patientID"].map(lambda x: "grp1" if x in grp1_patients else "grp2")
            filtered_data.obs = filtered_data.obs.rename(columns={
            "disease": "old_disease",
            "random_group": "disease"})
            filtered_data.obs["disease"] = filtered_data.obs["disease"].replace({"grp1": "CTL", "grp2": "AD"})

            # Save shuffled data
            output_path = os.path.join(output_dir, filename)
            filtered_data.write_h5ad(output_path, compression="gzip")
            print(f"Saved shuffled data to: {output_path}\n")

            uniq = filtered_data.obs.drop_duplicates(subset='patientID')

            # Now, create summary table
            summary = {
                "study_name": filename.split('_')[1].split('.')[0],
                "cell_type": filename.split('_')[0],
                "category": categories,
                "true_AD": np.sum(uniq["old_disease"] == "AD"),
                "true_CTL": np.sum(uniq["old_disease"] == "CTL"),
                "new_AD": np.sum(uniq["disease"] == "AD"),
                "new_CTL": np.sum(uniq["disease"] == "CTL"),
                "true_AD_in_new_AD": np.sum((uniq["old_disease"] == "AD") & (uniq["disease"] == "AD")),
                "true_CTL_in_new_AD": np.sum((uniq["old_disease"] == "CTL") & (uniq["disease"] == "AD")),
                "true_AD_in_new_CTL": np.sum((uniq["old_disease"] == "AD") & (uniq["disease"] == "CTL")),
                "true_CTL_in_new_CTL": np.sum((uniq["old_disease"] == "CTL") & (uniq["disease"] == "CTL"))
            }

            summary_df = pd.DataFrame([summary])
            summary_file = os.path.join(output_dir, f"{filename.split('.')[0]}_summary.csv")
            summary_df.to_csv(summary_file, mode='a', header=not os.path.exists(summary_file), index=False)

if __name__ == "__main__":
    main()


# python shuffle_patients_tr.py 
# --input_dir /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease 
# --output_dir /cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/01shuffledDisease/AD-Ctl
# --categories AD,CTL
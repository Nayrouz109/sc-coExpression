

import os
import pandas as pd


data_dir = "/cosmos/data/project-data/NW-rewiring/data/2.prcsd-slcGenes-cpm/intrsct-lessStringent-diag-disease/01shuffledDisease/AD-Ctl/tr"
summary_files = [f for f in os.listdir(data_dir) if f.endswith("_summary.csv")]
dfs = []

# Read each summary file and append it to the list
for file in summary_files:
    file_path = os.path.join(data_dir, file)
    df = pd.read_csv(file_path)
    df["Source_File"] = file  # Add a column to track the source file
    dfs.append(df)

# Concatenate all dataframes into one
if dfs:
    combined_df = pd.concat(dfs, ignore_index=True)
    
    # Save the concatenated dataframe to a CSV file
    output_file = os.path.join(data_dir, "combined_summary.csv")
    combined_df.to_csv(output_file, index=False)
    print(f"Combined summary saved to: {output_file}")
else:
    print("No summary files found in the directory.")

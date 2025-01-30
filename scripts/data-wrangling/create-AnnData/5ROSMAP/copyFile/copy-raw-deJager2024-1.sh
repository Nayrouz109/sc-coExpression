#!/bin/bash


# Define source and destination directories
src_dir="/cosmos/data/downloaded-data/AD-meta/ROSMAP-DeJager/data/0-source/expr/rds"
dest_dir="/space/scratch/nairuz-rewiring/data/0.processed-raw"

# List of files to copy (excluding the Exc file)
files=("Ast_DeJager.rds" \
       "Exc1_DeJager.rds" \
       "Exc2_DeJager.rds" \
       "Inh_DeJager.rds" \
       "Mic_DeJager.rds" \
       "Oli_DeJager.rds")

# Copy and rename files
for file in "${files[@]}"; do
    # Extract the cell type prefix (first part of the filename)
    cell_type=$(echo "$file" | cut -d'_' -f1)

    # Create the new filename with the required pattern
    new_filename="${cell_type}_ROSMAP-DeJager-2024.rds"

    # Copy the file to the destination and rename it
    cp "$src_dir/$file" "$dest_dir/$new_filename"
    
    # Print a message indicating success
    echo "Copied and renamed: $file to $new_filename"
done

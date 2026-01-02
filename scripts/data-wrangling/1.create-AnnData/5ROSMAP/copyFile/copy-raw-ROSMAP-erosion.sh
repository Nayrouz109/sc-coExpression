
#!/bin/bash

# Define source and destination directories
src_dir="/cosmos/data/downloaded-data/AD-meta/ROSMAP-erosion/data/0-source/expr/cell-type-rds"
dest_dir="/space/scratch/nairuz-rewiring/data/0.processed-raw"

# List of files to copy (excluding the Exc file)
files=("Ast_ROSMAP-erosion.rds" \
       "Exc_ROSMAP-erosion.rds" \
       "Inh_ROSMAP-erosion.rds" \
       "Mic_ROSMAP-erosion.rds" \
       "Opc_ROSMAP-erosion.rds" \
       "Oli_ROSMAP-erosion.rds")

# Copy and rename files
for file in "${files[@]}"; do
    # Extract the cell type prefix (first part of the filename)
    cell_type=$(echo "$file" | cut -d'_' -f1)

    # Create the new filename with the required pattern
    new_filename="${cell_type}_ROSMAP-erosion-2023.rds"

    # Copy the file to the destination and rename it
    cp "$src_dir/$file" "$dest_dir/$new_filename"
    
    # Print a message indicating success
    echo "Copied and renamed: $file to $new_filename"
done


#!/bin/bash

### Kamath-NPH ============================================================
# Define source and destination directories
src_dir="/cosmos/data/downloaded-data/AD-meta/Kamath-NPH/data/0-source/expr"
dest_dir="/space/scratch/nairuz-rewiring/data/0.processed-raw"

# List of files to copy (excluding the Exc file)
files=("Ast_data_arranged_updatedId_final_batches.rds" \
       "Exc_data_arranged_updatedId_final_batches.rds" \
       "Inh_data_arranged_updatedId_final_batches.rds" \
       "Mic_data_arranged_updatedId_final_batches.rds" \
       "Oli_data_arranged_updatedId_final_batches.rds" \
       "Opc_data_arranged_updatedId_final_batches.rds")

# Copy and rename files
for file in "${files[@]}"; do
    # Extract the cell type prefix (first part of the filename)
    cell_type=$(echo "$file" | cut -d'_' -f1)

    # Create the new filename with the required pattern
    new_filename="${cell_type}_Kamath-NPH-2023.rds"

    # Copy the file to the destination and rename it
    cp "$src_dir/$file" "$dest_dir/$new_filename"
    
    # Print a message indicating success
    echo "Copied and renamed: $file to $new_filename"
done

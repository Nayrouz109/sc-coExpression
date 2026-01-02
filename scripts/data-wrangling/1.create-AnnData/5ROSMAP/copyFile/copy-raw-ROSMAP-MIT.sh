
#!/bin/bash

# Define source and destination directories
src_dir="/cosmos/data/downloaded-data/AD-meta/ROSMAP-MIT/data/0-source/expr"
dest_dir="/space/scratch/nairuz-rewiring/data/0.processed-raw"

# List of files to copy
files=("Astrocytes.rds" \
       "Excitatory_neurons_set1.rds" \
       "Excitatory_neurons_set2.rds" \
       "Excitatory_neurons_set3.rds" \
       "Inhibitory_neurons.rds" \
       "Microglia_cells.rds"\
       "Oligodendrocytes.rds"\
       "OPCs.rds")

# Copy and rename files
for file in "${files[@]}"; do
    # Extract the cell type prefix (first part of the filename)
    cell_type=$(echo "$file" | cut -d'_' -f1)

    # Handle specific renaming based on filename
    if [[ "$file" == Excitatory_neurons_set* ]]; then
        # Extract set number for excitatory neurons
        set_num=$(echo "$file" | grep -o 'set[0-9]')
        new_filename="Exc${set_num: -1}_ROSMAP-MIT-2023.rds"
    else
        # Take first three letters of cell type for other files
        short_cell_type=$(echo "$cell_type" | cut -c1-3)
        new_filename="${short_cell_type}_ROSMAP-MIT-2023.rds"
    fi

    # Copy the file to the destination and rename it
    cp "$src_dir/$file" "$dest_dir/$new_filename"
    
    # Print a message indicating success
    echo "Copied and renamed: $file to $new_filename"
done

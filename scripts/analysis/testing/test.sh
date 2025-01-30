#!/usr/bin/bash

#SBATCH --job-name=test_hello       # Job name
#SBATCH --output=logs/test_hello_%j.out  # Standard output log
#SBATCH --error=logs/test_hello_%j.err   # Error log
#SBATCH --ntasks=1                  # Run one task
#SBATCH --cpus-per-task=1           # Use one CPU
#SBATCH --mem=1G                    # Allocate 1 GB of memory
#SBATCH --time=00:05:00             # Set a 5-minute time limit


# Activate the Conda environment
source ~/anaconda3/bin/activate homl3

# Print "Hello" using Python
python -c "print('Hello from the homl3 environment!')"

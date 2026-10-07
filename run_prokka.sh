#!/bin/bash

# Step 1: Create Conda environment with Prokka
echo "Creating Conda environment 'prokka'..."
conda create -y -n prokka -c conda-forge -c bioconda prokka

# Step 2: Activate environment
echo "Activating Conda environment 'prokka'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate prokka

# Step 3: Run Prokka on all .fna files in the current directory
echo "Running Prokka on all .fna files in $(pwd)..."
current_dir=$(pwd)

for file in "$current_dir"/*.fna; do
    filename=$(basename "$file" .fna)
    echo "Processing $filename..."
    prokka --outdir "$current_dir/$filename" --cpus 12 --prefix "$filename" "$file"
    echo "Prokka finished for $filename"
done

echo "All jobs completed."


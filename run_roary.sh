#!/bin/bash

# Step 1: Create Conda environment with Roary
echo "Creating Conda environment 'roary_env'..."
conda create -y -n roary_env -c conda-forge -c bioconda roary

# Step 2: Activate environment
echo "Activating Conda environment 'roary_env'..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate roary_env

# Step 3: Run Roary
echo "Running Roary on all .gff files..."
roary -e --mafft -p 12 *.gff

# Step 4: Download roary_plots.py script
echo "Downloading roary_plots.py..."
wget -q https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py

# Step 5: Install Python plotting libraries
echo "Installing matplotlib, seaborn, biopython..."
conda install -y matplotlib seaborn biopython

# Step 6: Run roary_plots.py (assumes output files exist)
echo "Running roary_plots.py..."
python roary_plots.py accessory_binary_genes.fa.newick gene_presence_absence.csv

echo "Roary analysis and plotting complete."


#!/bin/bash

# Create output folder if it doesn't exist
output_dir="gff_output"
mkdir -p "$output_dir"

# Find and copy all .gff files from subdirectories to gff_output
echo "Collecting .gff files..."
find . -type f -name "*.gff" -exec cp -u {} "$output_dir" \;

echo "All .gff files have been collected into '$output_dir'."


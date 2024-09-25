#!/bin/bash

# Loop through all .gb files in the test_species directory
for file in test_species/*.gb; do
  # Extract the base filename without the path and extension
  base=$(basename "$file" .gb)
  
  # Create the output directory based on the base name
  output_dir="results/$base"
  
  # Run the epsSMASH command with the file and output directory
  epsSMASH "$file" --output-dir "$output_dir"
  
  # Print the status
  echo "Processed $file, results saved to $output_dir"
done

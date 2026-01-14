#!/bin/bash

# Define the source and ground truth directories
source_dir="./build/unittest"
ground_truth_dir="./benchmarks/hrg/two_haplo/ground_truth/"

# Declare an associative array to track processed combinations
declare -A processed_combinations

# Iterate over all files matching the pattern in the source directory
for output_file in "$source_dir"/output_*.sam; do
  # Check if the file exists to avoid errors in case no files match
  if [[ -e "$output_file" ]]; then
    # Extract the filename from the full path
    filename=$(basename "$output_file")

    # Extract components from the output filename
    IFS='_' read -r testtype scheme part metric K I sparse i mode <<< "${filename#output_}"
    mode="${mode%.sam}" # Remove the .sam extension

    # Create a unique key for the combination of pair, name, part, metric, and I and mode!
    combination_key="${testtype}_${scheme}_${part}_${metric}_${I}_${mode}"

    # Check if this combination has already been processed
    if [[ -z "${processed_combinations[$combination_key]}" ]]; then
      # Mark this combination as processed
      processed_combinations[$combination_key]=1

      # Construct the new filename and path for ground truth
      ground_truth_file="${ground_truth_dir}${metric}${I}.${mode}.${testtype}.${scheme}.sam"

      # Move the output file to the ground truth directory
      mv "$output_file" "$ground_truth_file"

      echo "Moved $output_file to $ground_truth_file"
      echo "Combination: $combination_key"
    else
      echo "Skipped $output_file as it has the same combination as an already processed file."
      echo "Combination: $combination_key"
    fi
  else
    echo "No files found matching the pattern in $source_dir."
  fi
done

#!/bin/bash

# Check if the input directory is provided as a command-line argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

input_directory="$1"

# Extract the name of the input directory
input_directory_name=$(basename "${input_directory}")

# Construct the output directory based on the input directory name and tag
output_directory="${input_directory_name}_root"

# Create the output directory if it doesn't exist
mkdir -p "${output_directory}"

# Specify the name for the final merged root file based on the input directory name
merged_root_file="${input_directory_name}_merged.root"

# Get the list of LHE files in the input directory
lhe_files=("${input_directory}"/pwgevents-*.lhe)

# Loop through all the LHE files
for lhe_file_path in "${lhe_files[@]}"; do
    # Extract the file number from the LHE file name
    file_number=$(basename "${lhe_file_path}" | sed -n 's/^pwgevents-\([0-9]\+\)\.lhe$/\1/p')

    # Construct the output root file name
    root_file="pwgevents-${file_number}.root"

    # Remove the line containing the specified pattern because it breaks the reading script
    sed -i '/#Random number generator exit values:/d' "${lhe_file_path}"

    # Run your Python script to convert LHE to root
    python3 LHEReader.py --input "${lhe_file_path}" --output "${output_directory}/${root_file}"

    echo "Processed: ${lhe_file_path} -> ${output_directory}/${root_file}"
done

# Merge all individual root files into one
hadd "${merged_root_file}" "${output_directory}"/*.root 

echo "All root files are merged into ${merged_root_file}"

# Specify the directory for the plots
plots_directory="${input_directory_name}_plots"

# Run your Python script to create plots from the merged root file
python3 LHEplots.py --input "${merged_root_file}" --output "${plots_directory}"

echo "Plots created in ${plots_directory}"
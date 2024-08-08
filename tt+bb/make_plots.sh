#!/bin/bash

if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_root_file>"
    exit 1
fi

input_root_file="$1"

# Check if the input_root_file ends with "_merged.root"
if [[ $input_root_file == *_merged.root ]]; then
    # Remove "_merged.root" and store the base name
    base_name="${input_root_file%_merged.root}"
else
    echo "Error: Input file does not end with '_merged.root'"
    exit 1
fi

output_directory="${base_name}_plots"

python3 LHEplots_FH.py --input "${input_root_file}" --output "${output_directory}"
# python3 LHEplots_SL.py --input "${input_root_file}" --output "${output_directory}"
# python3 LHEplots_DL.py --input "${input_root_file}" --output "${output_directory}"

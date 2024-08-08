#!/bin/bash

# A bash script that runs the ttbb.C script for multiple files
# In order to run it, change the number in the for loop to the amount of the root files you have and of course the name of the files.

for i in {1..10}; do
    input_file="lhe_hdampUP_FH_${i}__r1.0_f1.0_m172.5_p320900_merged.root"
    output_file="lhe_hdampUP_FH_${i}__r1.0_f1.0_m172.5_p320900_ttbb_plots.root"
    root -l -b -q "ttbb.C(\"$input_file\", \"$output_file\")"
done


# Scripts for tt+bb 
Scripts that read Les Houches Event (LHE) files and create plots for the tt+bb process

## Usage

In order to use the scripts for the tt+bb process simply replace the `LHEplots.py` script with the desired one. Assuming that the `.lhe` has already been converted into a `.root` file, run:

```
python3 LHEplots_[CHANNEL].py --input [ROOT_FILE] --output [FOLDER]
```

**NOTE**: If you have a very large file it is recommended to use a C++ script instead of the python ones. 

An example that calculates the invariant mass and the HT of the tt+bb system is provided (same for all three channels). In order to execute it, simply run:

```
root -l -b -q "ttbb.C(\"[INPUT_FILE]\", \"[OUTPUT_FILE]\")"
```

If you want to run it for multiple files, use the corresponding bash script:
```
./ttbb.sh
```

## Overlay 

The above scripts, other than the plots(in both `.png` and `.pdf` format), create also a ROOT file that contains all the histograms in a tree. Using this file we can further process it and overlay it with others in order to compare different settings and etc. For this purpose the `ratio_[CHANNEL].py` scripts have been created. In order to use them, simply run:
```
python3 ratio_[CHANNEL].py
```
and a folder will be created with all the plots.  

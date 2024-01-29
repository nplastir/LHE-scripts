# LHE-scripts
Python scripts that read Les Houches Event (LHE) files and create plots 

## Setup

```
git clone https://github.com/nplastir/LHE-scripts.git
```

## Usage 

- Step 1: Read LHE file and convert it in a ROOT tree

```
python3 LHEReader.py --input input.lhe --output output.root
```
- Step 2: Read ROOT tree and create plots
```
python3 LHEplots.py --input output.root --output [DIRECTORY]
```


## Example
```
python3 LHEReader.py --input example.lhe --output example.root
```
Running this, will create a root file ```example.root```. Then, simply run
```
python3 LHEplots.py --input example.root --output example
```
which will create a directory (if does not already exist) called ```example``` where will be created some plots in ```.png``` format. 


## Batch conversion

If you have multiple files that you want to convert into ROOT files:
- Create a folder where you will transfer all the ```.lhe``` files that you want to process 
- Specify the path to the directory containing your LHE files that you have created inside the ```batch_convert.sh``` script and run it ```./batch_convert.sh```
- This script will then create, according to the input folder name, a directory to store all the ROOT files, will merge them and then create a directory with all the plots 

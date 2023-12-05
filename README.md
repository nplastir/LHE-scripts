# LHE-scripts
Python scripts that read Les Houches Event (LHE) files and create plots 

## Setup

```
git clone https://github.com/nplastir/LHE-scripts.git
```

## Usage 

- Step 1: Read LHE file and convert it in a ROOT tree

```
python3 LHEReader --input input.lhe --output output.root
```
- Step 2: Read ROOT tree and create plots:
```
python3 LHEplots.py --input output.root --output [DIRECTORY]
```


## Example: 
```
python3 LHEReader.py --input example.lhe --output example.root
```
Running this, will create a root file ```example.root```. Then, simply run
```
python3 LHEplots.py --input example.root --output example
```
which will create a directory (if does not already exist) called ```example``` where will be created some plots in ```.png``` format. 
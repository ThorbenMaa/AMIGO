# AMIGO
Here you will find information on how to set up AMIGO on a Linux or Mac (Sonoma 14.1.1) system. Make sure you have installed git. If not, do so by using the follwoing command:
```
sudo apt-get install git
```

### Clone this repository
Use 
```
git clone https://github.com/ThorbenMaa/AMIGO.git
```
Operate from inside the directory.

### Install dependencies
I recommend to use [mamba](https://mamba.readthedocs.io/en/latest/installation.html) to create environments and install dependencies.
After setting up mamba or, alternatively, conda, set up the environment using:

```
mamba env create --name AMIGO --file=./env/AMIGO.yml
mamba activate AMIGO
```
or
```
mamba create -n AMIGO
mamba activate AMIGO
mamba install python
mamba install networkx
mamba install numpy
mamba install matplotlib
mamba install pandas
mamba install click
```

### Run AMIGO
For running AMIGO, please look at the files in `InputFilesExample_LmUGP` or in `InputFilesExample_GTB`. Further descriptions of the corresponding files can be found in Supplementary Material of the manuscript. Submitting parameters to AMIGO can be done in two ways: 
- interactively by placing and naming the files `noe.txt`, `preassignments.txt`, `additional_measure.txt`, `additional_xtal.txt`, and `AMIGO.py` in one folder as in the example folders and using the command `python AMIGO.py` from inside this folder.
- placing `run_AMIGO.sh` and `AMIGO_non_interactive.py` in one folder, adjusting the parameters and specifying files with their respective paths in `run_AMIGO.sh`, and using the command `bash run_AMIGO.sh` from inside this folder.


By running `bash run_AMIGO.sh`, `AMIGO_non_interactive.py` is executed. By running `python AMIGO.py`, the interactive version of AMIGO is directly executed. Both are an identical implementation of the methyl walk steps described in our manuscript (REF). For more information on input files, on optimization and choice of input parameters, and on output files, see the SI information of our manuscript.


Note that the script `AMIGO_corr.py` in subfolders of the `Benchmark` folder uses the same implementation as in `AMIGO.py` or `AMIGO_non_interactive.py` in the example folders. The only difference is that there are some typos/changes in the documentation and printing of the result files. We provide `AMIGO_corr.py` for the sake of completeness and to be able to exactly reproduce the files in the `Benchmark` folder. Parameters can be submitted interactively using the ones specified in the respective results.txt files.

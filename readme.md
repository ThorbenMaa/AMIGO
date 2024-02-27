
# install AMIGO
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

# Tutorial on how to use AMIGO without PCSs on the LmUGP example on a Linux system
## cut-off range optimization
### prepare cut-off range optimization
Operate from inside the `Tutorial_no_PCSs_LmUGP` folder.

If you have your own data set, make sure that the data folder contains the according `noe.txt`, `additional_measure.txt`, `additional_xtal.txt`, and `preassignments.txt` file as in the `data` folder. `additional_measure.txt`, `additional_xtal.txt`, and `preassignments.txt` are optional and you can leave them empty. However, the (empty) files have to be present. We won't use any information in `additional_measure.txt` or `additional_xtal.txt` in this tutorial by setting the according weight factors to 0 (see below).

**Details of the noe.txt file**
The first column contains the amino acid type of the NOE donor, the second column contains the peak ID of the NOE donor,
the third column contains the amino acid type of the NOE acceptor, the fourth column contains the peak ID of the NOE acceptor. **The file must be sorted using the second column**. If you have no stereospecific Leu and Val labeling, and there is more than one NOE between a Leu or Val residue and another particular residue (i.e., the proR and the proS methyl group have an NOE with to the same neighbor), use this NOE only once. E.g.:
```
23 LEU_proS 29 ALA
23 LEU_proR 29 ALA
```
would become

```
23 LEU 29 ALA
```
or 
```
29 ALA 23 LEU_proS
29 ALA 23 LEU_proR
```
would become
```
29 ALA 23 LEU
```

Generate a copy of the `data` folder for each cut-off range you want to test. The test runs will be run with a reduced number of possible methyl walk starting points to make the calculations quick and should not be used for generating a final assignment.
```
cp -R data 7
cp -R data 7-5
cp -R data 8
cp -R data 8-5
cp -R data 9
cp -R data 9-5
cp -R data 10
```

Now adopt the `run_AMIGO.sh` files within the generated folders according to the example for `7` given here:
```
#!/bin/bash

python AMIGO_non_interactive.py \
--pdb_file 4m2a_pro.pdb \
--noe_file noe.txt \
--additional_measure_file additional_measure.txt \
--additional_xtal_file additional_xtal.txt \
--preassignments_file preassignments.txt \
--min_met 3 \
--max_met 7 \
--step_met 0.2 \
--min_ile 3 \
--max_ile 7 \
--step_ile 0.2 \
--min_leu 3 \
--max_leu 7 \
--step_leu 0.2 \
--min_val 3 \
--max_val 7 \
--step_val 0.2 \
--min_ala 3 \
--max_ala 7 \
--step_ala 0.2 \
--min_thr 3 \
--max_thr 7 \
--step_thr 0.2 \
--weight_noe 1 \
--weight_additional_resraint1 0 \
--weight_additional_resraint2 0 \
--weight_additional_resraint3 0 \
--weight_additional_resraint4 0 \
--percentage_starting_points 0.1 \
--val_scheme proS \
--leu_scheme proS
```

For the other folders (`7-5`, `8`, `8-5`, ...) change the maximal cut-off distance parameters accordingly, too (e.g., for `10` use `max_met 10`, `max_ile 10`, `max_leu 10`, ...). Make also sure to set `percentage_starting_points` to `0.1` and to choose the right labeling scheme for val and leu. If you have racemic mixtures use `both`. If you have proR labeling use `proR`. If you haven't labeled a particular amino acid type, e.g., THR, use `0`, `0`, `1` for `min_thr`, `max_thr`, and `step_thr`, respectively.


### run cut-off range optimization
To run the cut-off optimization, run
```
cd ./7
bash AMIGO_non_interactive.sh
cd ./../7-5
bash AMIGO_non_interactive.sh
cd ./../8
bash AMIGO_non_interactive.sh
cd ./../8-5
bash AMIGO_non_interactive.sh
cd ./../9
bash AMIGO_non_interactive.sh
cd ./../9-5
bash AMIGO_non_interactive.sh
cd ./../10
bash AMIGO_non_interactive.sh

```

### determine optimal cut-off ranges
If all jobs are finished you can extract the number of total assignments and the number of perfectly matching building blocks by using

```
# total number of assignments
tail ./7/result.txt | grep "total amount of assignments"
# number of perfectly matching building blocks
tail ./7/result.txt | grep "total amount of perfectly matching building block pairs"
```

and calculate the ratio (e.g. in an excel or libre office file). Replace the "7" with "7-5", "8", "8-5", ... to get the ratios for all cut-off ranges.


The three best cut-off ranges based on the ratio should be the ones having 7, 7.5, and 8 Angstroem as maximal cut-off distance. The best one should be 7.5 Angstroem.

## actual AMIGO runs
### prepare the folders for the actual runs
Create a folder for the actual AMIGO runs.
```
mkdir p=1
```

and copy the folders corresponding to the best cut-off ranges into the new folder:

```
cp -R 7 p\=1/7
cp -R 7-5 p\=1/7
cp -R 8 p\=1/7
```

Now adopt the `run_AMIGO.sh` files within the generated sub folders in `p=1` according to the example for `7` given here:
```
#!/bin/bash

python AMIGO_non_interactive.py \
--pdb_file 4m2a_pro.pdb \
--noe_file noe.txt \
--additional_measure_file additional_measure.txt \
--additional_xtal_file additional_xtal.txt \
--preassignments_file preassignments.txt \
--min_met 3 \
--max_met 7 \
--step_met 0.2 \
--min_ile 3 \
--max_ile 7 \
--step_ile 0.2 \
--min_leu 3 \
--max_leu 7 \
--step_leu 0.2 \
--min_val 3 \
--max_val 7 \
--step_val 0.2 \
--min_ala 3 \
--max_ala 7 \
--step_ala 0.2 \
--min_thr 3 \
--max_thr 7 \
--step_thr 0.2 \
--weight_noe 1 \
--weight_additional_resraint1 0 \
--weight_additional_resraint2 0 \
--weight_additional_resraint3 0 \
--weight_additional_resraint4 0 \
--percentage_starting_points 1 \
--val_scheme proS \
--leu_scheme proS
```

The only parameter that is different is `percentage_starting_points` which was `0.1` before and which we changed to `1`. Do the same thing for `7-5` and `8` accordingly.

### start the actual runs
Strat the actual runs by navigating into the `p=1` folder  
```
cd p\=1
```

and using

```
cd ./7
bash AMIGO_non_interactive.sh
cd ./../7-5
bash AMIGO_non_interactive.sh
cd ./../8
```

### analyse the results
If you want to go through the methyl walks yourself, it is recommended to use the run with the optimal cut-off range (which would be 7-5 based on the ratio calculated in the cut-off range optimization steps). If you rather prefer only very safe assignments, you should only use assignments shared in the runs with the three top cut-off ranges (see above). To do so, you can use the `compare.py` script in the Tutorial folder.

```
# copy compare script
cp ./../compare.py p\=1

# extract shared assignments
python compare.py \
--results_files 7/results.txt \
--results_files 7-5/results.txt \
--results_files 8/results.txt \
--outfile_name shared_assignments.tsv
```

The shared assignments will end up in the file `shared_assignments.tsv`.

The folder structure and files should look identical to `LmUGP_without_PCS` in the `Benchmark/benchmark_data_set` folder.

# Tutorial on how to use AMIGO with PCSs on the GTB example on a Linux system
Operate from inside the `Tutorial_PCS_GTB` folder.

Using PCS makes AMIGO very robust in terms of cut-off range choices. We recomment to just use 3 to 9 Angstroems in 0.2 Angstroem steps for all labeled amino acids. Prepare all input files according to the files in this folder.
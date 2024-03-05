- [install AMIGO](#install-amigo)
    - [Clone this repository](#clone-this-repository)
    - [Install dependencies](#install-dependencies)
- [Tutorial on how to use AMIGO without PCSs on the LmUGP example on a Linux system](#tutorial-on-how-to-use-amigo-without-pcss-on-the-lmugp-example-on-a-linux-system)
  - [cut-off range optimization](#cut-off-range-optimization)
    - [prepare cut-off range optimization](#prepare-cut-off-range-optimization)
    - [run cut-off range optimization](#run-cut-off-range-optimization)
    - [determine optimal cut-off ranges](#determine-optimal-cut-off-ranges)
  - [actual AMIGO runs](#actual-amigo-runs)
    - [prepare the folders for the actual runs](#prepare-the-folders-for-the-actual-runs)
    - [start the actual runs](#start-the-actual-runs)
    - [analyse the results](#analyse-the-results)
- [Tutorial on how to use AMIGO with PCSs on the GTB example on a Linux system](#tutorial-on-how-to-use-amigo-with-pcss-on-the-gtb-example-on-a-linux-system)
  - [prepare the files](#prepare-the-files)
    - [start the actual runs](#start-the-actual-runs-1)
    - [analyse the results](#analyse-the-results-1)

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

If you have your own data set, make sure that the data folder contains the according `noe.txt`, `additional_measure.txt`, `additional_xtal.txt`, and `preassignments.txt` file as in the `data` folder. `additional_measure.txt`, `additional_xtal.txt`, and `preassignments.txt` are optional. However, the files have to be present. For `additional_measure.txt`, `additional_xtal.txt` you can put a single line in there containing `0	0	0	0	0	0`. `preassignments.txt` can be left completely empty.

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
If you want to go through the methyl walks yourself, it is recommended to use the run with the optimal cut-off range (here, these methyl walks would be in the `result.txt` file in the `7-5` folder (based on the ratio calculated in the cut-off range optimization steps)). The files `G_NMR.svg` and `G_pdb.svg` contain the NOE and structure-based networks. You can try to follow the methyl walks in `result.txt` in these networks.


If you want to consider only very safe assignments, you should only use assignments shared in the runs with the three top cut-off ranges (see above). To do so, you can use the `compare.py` script in the Tutorial folder.

```
# copy compare script
cp ./../compare.py p\=1

# extract shared and non-shared assignments and unassigned peak IDs
python compare.py \
--results_files 7/results.txt \
--results_files 7-5/results.txt \
--results_files 8/results.txt \
--min_cluster_size 1 \
--noe_file 7/noe.txt \
--outfile_name shared_assignments.tsv
```

The shared assignments will end up in the file `shared_assignments.tsv`. Ambiguous assignments not assigned to the same amino acid in the three runs will end up in `ambiguous_shared_assignments.tsv`. Not assigned peak IDs will end up in `unassigned_shared_assignments.tsv`.

The folder structure and files should look similar to `LmUGP_without_PCS` in the `Benchmark/benchmark_data_set` folder.

# Tutorial on how to use AMIGO with PCSs on the GTB example on a Linux system
Operate from inside the `Tutorial_PCS_GTB` folder.


## prepare the files
Using PCS makes AMIGO very robust in terms of cut-off range choices. We recomment to use three runs using 3 Angstroem as minimal cut-off distance, 8, 8.5 and 9 Angstroem as maximal cut-off distance, and a step size of 0.2 Angstroem for all labeled amino acids. Prepare all input files according to the files in the `data` folder. 

If you have your own data set, make sure that the data folder contains the according `noe.txt` (explanations given in the tutorial on how to use AMIGO without PCSs), `additional_measure.txt`, `additional_xtal.txt`, and `preassignments.txt` file as in the `Tutorial_PCS_GTB` folder. `preassignments.txt` are optional and you can leave it empty. However, the (empty) file has to be present. 

**Details of the `additional_measure.txt`, `additional_xtal.txt`**
`additional_measure.txt` contains measured PCSs, `additional_xtal.txt` contains theoretically determined PCSs. Up to four different sets of PCSs are supported (e.g. induced by different paramagnetic metals). The first and second columns are the amino acid type and the peak ID or the amino acid type and the residue ID from the pdb file, respectively. The third columns of the two files correspond to the first set of measured and theoretically derived PCSs, respectively. The second ones to the second set of measured and theoretically derived PCSs, respectively, and so on... If you have no values for a given set (e.g., the last column in the example files) use `999` as a placeholder. If you haven't measured a value for a particular methyl group, e.g., because it is close to the paramagnetic center and the peak is broadened beyond detection, also use `999` (see also example file, e.g.,  MET 69).

Make three copies of the data folder for the three different runs:

```
cp -R data 8
cp -R data 8-5
cp -R data 9
```

Now adopt the `run_AMIGO.sh` files within the generated folders according to the example for `8` given here:
```
#!/bin/bash

python AMIGO_non_interactive.py \
--pdb_file 4m2a_pro.pdb \
--noe_file noe.txt \
--additional_measure_file additional_measure.txt \
--additional_xtal_file additional_xtal.txt \
--preassignments_file preassignments.txt \
--min_met 3 \
--max_met 8 \
--step_met 0.2 \
--min_ile 3 \
--max_ile 8 \
--step_ile 0.2 \
--min_leu 3 \
--max_leu 8 \
--step_leu 0.2 \
--min_val 3 \
--max_val 8 \
--step_val 0.2 \
--min_ala 3 \
--max_ala 8 \
--step_ala 0.2 \
--min_thr 0 \
--max_thr 0 \
--step_thr 1 \
--weight_noe 1 \
--weight_additional_resraint1 100 \
--weight_additional_resraint2 100 \
--weight_additional_resraint3 100 \
--weight_additional_resraint4 0 \
--percentage_starting_points 0.1 \
--val_scheme proS \
--leu_scheme proS
```

For the other folders (`8-5` and `9`) change the maximal cut-off distance parameters accordingly, too (e.g., for `9` use `max_met 9`, `max_ile 9`, `max_leu 9`, ...). Make also sure to set `percentage_starting_points` to `1` and to choose the right labeling scheme for val and leu. If you have racemic mixtures use `both`. If you have proR labeling use `proR`. If you haven't labeled a particular amino acid type, e.g., THR, use `0`, `0`, `1` for `min_thr`, `max_thr`, and `step_thr`, respectively. Also note that we have set the weight factors for additional restraints `weight_additional_resraint1`, `weight_additional_resraint2`, and `weight_additional_resraint3` to `100`. `100` is a good value for PCSs used here. The differences between measured and theoretical PCSs times the according weight factor should be in the 10^0 to 10^1 range.

### start the actual runs
Strat the actual runs by using:

```
cd ./8
bash AMIGO_non_interactive.sh
cd ./../8-5
bash AMIGO_non_interactive.sh
cd ./../9
```

### analyse the results
You can find the methyl walks annotated with the corresponding PCSs in the `result.txt` files of the respecitve folders (`8`, `8-5`, and `9`). The folders will also contain a file named `additional_restraints.svg` showing correlation plots of theoretical and measured PCSs based on the assignment. If the assignment is reasonably good, they should be linaerly correlated. The files `G_NMR.svg` and `G_pdb.svg` contain the NOE and structure-based networks. You can try to follow the methyl walks in `result.txt` in these networks. If you want to consider only very safe assignments, you should only use assignments shared in the three runs. To do so, you can use the `compare.py` script in the Tutorial folder.

```
# extract shared and non-shared assignments and unassigned peak IDs
python compare.py \
--results_files 8/results.txt \
--results_files 8-5/results.txt \
--results_files 9/results.txt \
--min_cluster_size 1 \
--noe_file 8/noe.txt \
--outfile_name shared_assignments.tsv
```

The shared assignments will end up in the file `shared_assignments.tsv`. Ambiguous assignments not assigned to the same amino acid in the three runs will end up in `ambiguous_shared_assignments.tsv`. Not assigned peak IDs will end up in `unassigned_shared_assignments.tsv`.

The folder structure and files should look similar to `GTB` in the `Benchmark/benchmark_data_set/PCS_data_sets/` folder (you can ignore the runs with 9.5 and 10 Angstroem as maximal cut-off distance in `9-5` and `10`).
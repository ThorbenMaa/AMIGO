for the noes generated here, the following command has been used with cut off disatnces 6.5, 7.5, and 8.5:

```
python simulate_NOE_network_without_random_removal.py --pdb_file ./../PCS_data_sets/LmUGP_with_PCS/8/Holo_run_1\(4m2a\).pdb \
--leu_scheme proS \
--val_scheme proS \
--dist_cut_off 6.5 \
--outfile ./noe_4m2a_6_5_A.txt
```
Note that using this script, NOEs are not randomly removed but the full networks are created. 
The outpput files were moved manually to this folder afterwards.
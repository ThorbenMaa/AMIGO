for the noes generated here, the following command has been used with cut off disatnces 6.5, 7.5, and 8.5 and with the closed an both open conformations:

```
python simulate_NOE_network_without_random_removal.py --pdb_file ./../PCS_data_sets/GTB/8/GTB_closed_conformation_final_model.pdb \
--leu_scheme proS \
--val_scheme proS \
--dist_cut_off 6.5 \
--outfile ./GBT_simulated_NOEs/noe_GTB_closed_conformation_final_model_6_5_A.txt

python simulate_NOE_network_without_random_removal.py \
--pdb_file ./../PCS_data_sets/GTB_conf_1_open_final_0_1_stepsize/8/GTB_conf_1_open_final.pdb \
--leu_scheme proS \
--val_scheme proS \
--dist_cut_off 6.5 \
--outfile ./GBT_simulated_NOEs/noe_GTB_conf_1_open_final_6_5_A.txt
number of methlys 194


```
Note that using this script, NOEs are not randomly removed but the full networks are created. 
The outpput files were moved manually to this folder afterwards.
#! /bin/bash
protein=$1
pdb=$2


python simulate_NOE_network.py \
--pdb_file ${protein}/${pdb} \
--leu_scheme proR \
--val_scheme proR \
--dist_cut_off 7 \
--outfile ${protein}/noe.txt




#!/bin/bash

python AMIGO_non_interactive.py \
--pdb_file 6y7f.pdb \
--noe_file noe.txt \
--additional_measure_file additional_measure.txt \
--additional_xtal_file additional_xtal.txt \
--preassignments_file preassignments.txt \
--min_met 3 \
--max_met 7.5 \
--step_met 0.2 \
--min_ile 3 \
--max_ile 7.5 \
--step_ile 0.2 \
--min_leu 3 \
--max_leu 7.5 \
--step_leu 0.2 \
--min_val 3 \
--max_val 7.5 \
--step_val 0.2 \
--min_ala 3 \
--max_ala 7.5 \
--step_ala 0.2 \
--min_thr 0 \
--max_thr 0 \
--step_thr 1 \
--weight_noe 1 \
--weight_additional_resraint1 0 \
--weight_additional_resraint2 0 \
--weight_additional_resraint3 0 \
--weight_additional_resraint4 0 \
--percentage_starting_points 1 \
--val_scheme proS \
--leu_scheme proS




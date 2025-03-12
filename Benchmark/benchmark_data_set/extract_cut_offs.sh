#!/bin/bash

# no PCS data
paths=("./HuNoV_Pdomain/p=1/9" \
        "./msg/p=1/8-5" \
        "./mbp/p=1/9" \
        "./atcase_fix/p=1/10" \
        "./EIN/p=1/7" \
        "./mrsB/p=1/9-5" \
        "./ubiquitin/p=1/7" \
        "./a7a7/p=1/10" \
        "./LmUGP_apo_without_PCS/p=1/7-5" \
        "./LmUGP_without_PCS/p=1/7-5" \
        "./GTB_without_PCS/p=1/8-5" \
        "./MNV_Pdomain_without_PCS/p=1/8" \
)

for path in "${paths[@]}"
do
    echo "Processing: ${path}"
    cat $path/result.txt | grep "and thereby assigned to pdb" | tr -s " " "\t" | awk '{print $9, $10, $(NF-1), $NF}' > ${path}/cut_off_dists.tsv
done

# PCS data
cd PCS_data_sets

paths=("GTB/8-5" "LmUGP_with_PCS/8-5" "MNV_Pdomain_with_PCS/8-5" "GTB_0_1_step_size_no_assignment_restrictions_only_loops/8-5" "GTB_assign_only_loops_0_1_step_size_no_PCSs/7-5")

for path in "${paths[@]}"
do
    echo "Processing: ${path}"
    cat $path/result.txt | grep "and thereby assigned to pdb" | tr -s " " "\t" | awk '{print $9, $10, $(NF-1), $NF}' > ${path}/cut_off_dists.tsv
done

cd ..
echo done
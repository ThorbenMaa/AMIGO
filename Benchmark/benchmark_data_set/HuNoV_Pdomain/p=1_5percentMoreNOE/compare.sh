#!/bin/bash
for i in 0 1 2 3 4 5 6;
do
  python compare.py \
  --noe_file 9/noe.txt \
  --results_files 7/result.txt \
  --results_files 8-5/result.txt \
  --results_files 9/result.txt \
  --min_cluster_size $i \
  --outfile_name shared_min_cluster_$i.tsv
done



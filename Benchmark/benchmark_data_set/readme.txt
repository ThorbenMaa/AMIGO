this folder contains subflders with the proteins used for benchamrking AMIGO. The results of cut-off distance range optimization, perfomance using the optimal cut-off range, and performance using the top three cut-off ranges are summarised in `analysis_final.ods`. 

Cut-off distance optimization of synthetic data sets can be found in `asuwertung_synthetic.ods`.

All data (except synthetic MILVATproR) sets including cut-off optimization can be run by using the `sbatch_run_benchmark_time.sh` script.

A performance comparison to other algorithms can be found in `comparison_to_other_algos_final.ods`.

Benchmarking results of data sets with PCS can be found in `PCS_data_sets`.

To run cut-off distance range optimiation for all data sets except of `HuNoV_Pdomain` and `MNV_Pdomain_without_PCS` and the datasets in `PCS_data_sets` you can use the script `run_benchmark_par_opt.sh`.

Also see readme files of the respective subfolders.

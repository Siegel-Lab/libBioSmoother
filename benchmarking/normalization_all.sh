#!/bin/bash

# while [[ $(squeue | grep $(whoami)) ]]
# do
#     echo "jobs running, waiting 60 seconds"
#     sleep 60
# done
# echo "no jobs running, starting benchmarking"

RADICL_INDEX=radicl-mouse/index.50000.smoother_index
HIC_INDEX=organisms/t_brucei/out/HGAP3_Tb427v10.10000.1.smoother_index

for bin_size in 50000 100000 500000
do
    window_size=0

    # grid seq ground truth
    out_file=norm_out/grid_seq_ground_truth.${bin_size}
    extra_param="settings.normalization.normalize_by=grid-seq settings.normalization.grid_seq_global=True"
    sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    # radicl seq ground truth
    out_file=norm_out/radicl_seq_ground_truth.${bin_size}
    extra_param="settings.normalization.normalize_by=radicl-seq settings.normalization.radicl_local=True"
    sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    # raw data (radicl)
    out_file=norm_out/raw_data_radicl.${bin_size}
    extra_param="settings.normalization.normalize_by=dont"
    sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    # ddd
    out_file=norm_out/ddd_ground_truth.${bin_size}
    extra_param="settings.normalization.normalize_by=dont settings.normalization.ddd=True settings.normalization.ddd_samples.val_min=0 settings.normalization.ddd_all_samples=True settings.filters.symmetry=mirror"
    sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    # cooler
    out_file=norm_out/cooler_ground_truth.${bin_size}
    extra_param="settings.normalization.normalize_by=cool-ice settings.filters.symmetry=mirror"
    sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    # local-ice
    out_file=norm_out/ice_ground_truth.${bin_size}
    extra_param="settings.normalization.ice_local=True settings.normalization.normalize_by=ice settings.filters.symmetry=mirror"
    sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    # raw data (hic)
    out_file=norm_out/raw_data_hic.${bin_size}
    extra_param="settings.normalization.normalize_by=dont settings.filters.symmetry=mirror"
    sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

    window_size=5000000
    for num_samples in 1 2 10 25 100 250 1000
    do
        # grid-seq
        out_file=norm_out/grid_seq.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=grid-seq settings.normalization.grid_seq_samples.val=${num_samples}"
        sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

        # radicl-seq
        out_file=norm_out/radicl_seq.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=radicl-seq settings.normalization.radicl_seq_samples.val=${num_samples}"
        sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

        # ice
        out_file=norm_out/ice.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=ice settings.normalization.num_ice_bins.val=${num_samples} settings.filters.symmetry=mirror"
        sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

        # ddd
        out_file=norm_out/ddd.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=dont settings.normalization.ddd=True settings.normalization.ddd_samples.val_min=0 settings.normalization.ddd_samples.val_max=${num_samples} settings.filters.symmetry=mirror"
        sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh
    done
done

bin_size=100000
for window_size in 1000000 5000000 10000000
do
    for num_samples in 1 2 10 25 100 250 1000
    do
        # grid-seq
        out_file=norm_out/grid_seq.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=grid-seq settings.normalization.grid_seq_samples.val=${num_samples}"
        sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

        # radicl-seq
        out_file=norm_out/radicl_seq.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=radicl-seq settings.normalization.radicl_seq_samples.val=${num_samples}"
        sbatch --export INDEX_FILE=${RADICL_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

        # ice
        out_file=norm_out/ice.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=ice settings.normalization.num_ice_bins.val=${num_samples} settings.filters.symmetry=mirror"
        sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh

        # ddd
        out_file=norm_out/ddd.${bin_size}.${window_size}.${num_samples}
        extra_param="settings.normalization.normalize_by=dont settings.normalization.ddd=True settings.normalization.ddd_samples.val_min=0 settings.normalization.ddd_samples.val_max=${num_samples} settings.filters.symmetry=mirror"
        sbatch --export INDEX_FILE=${HIC_INDEX},OUTPUT_FILE="${out_file}",BIN_SIZE=${bin_size},WINDOW_SIZE=${window_size},EXTRA_PARAMS="${extra_param}" bin/normalization.sh
    done
done
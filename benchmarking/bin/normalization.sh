#!/bin/bash

#SBATCH --cpus-per-task=18
#SBATCH --job-name=smoother_normalization
#SBATCH --mem=250G
#SBATCH --time=12-00:00:00
#SBATCH -o norm_out/slurm-out-normalization-%j.out
#SBATCH --partition=slim18

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


echo -n "Error or in progress" > ${OUTPUT_FILE}.status

echo ${INDEX_FILE} ${OUTPUT_FILE} ${BIN_SIZE} ${WINDOW_SIZE} ${EXTRA_PARAMS}
python -u bin/test_normalizations.py ${INDEX_FILE} ${OUTPUT_FILE} ${BIN_SIZE} ${WINDOW_SIZE} ${EXTRA_PARAMS}

echo -n "OK" > ${OUTPUT_FILE}.status
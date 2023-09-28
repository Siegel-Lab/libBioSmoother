#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=smoother_benchmark
#SBATCH --mem=250G
#SBATCH --partition=slim18
#SBATCH --time=12-00:00:00
#SBATCH -o out/slurm-out-smoother_benchmark-%j.out

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


GENOME_FILES=(genome/*.fna)
GENOME_FILE=$(echo ${GENOME_FILES[0]} | sed "s/.fna//")
INDEX_FILE=$(echo ${GENOME_FILE} | sed "s/genome/out/g")

EXTRA_SMOOTHER_PARAMS_FILE=$(echo ${EXTRA_SMOOTHER_PARAMS} | sed "s/ /./g")

FULL_SUFFIX=${RESO}.${SUBSAMPLE_FRACTION}${EXTRA_SMOOTHER_PARAMS_FILE}

echo "benchmarking index for ${GENOME_FILE} at resolution ${RESO}, subsamples ${SUBSAMPLE_FRACTION}, and extra_params ${EXTRA_SMOOTHER_PARAMS}"

INPUT_FOLDER=fasta
OUTPUT_FOLDER=data

QUERY_STATUS_FILE=out/index_query_times
QUERY_TIME_FILE=out/benchmark

for REPLICATE in {0..99}
do
    echo -n "Error or in progress" > ${QUERY_STATUS_FILE}.${FULL_SUFFIX}.${REPLICATE}.txt

    biosmoother benchmark ${INDEX_FILE}.${FULL_SUFFIX}.smoother_index -o ${QUERY_TIME_FILE}.${FULL_SUFFIX}.${REPLICATE}.pickle -N 1000

    echo -n "OK" > ${QUERY_STATUS_FILE}.${FULL_SUFFIX}.${REPLICATE}.txt
done
echo "done"
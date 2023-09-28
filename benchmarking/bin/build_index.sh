#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=smoother_index
#SBATCH --mem=250G
#SBATCH --partition=slim18
#SBATCH --time=12-00:00:00
#SBATCH -o out/slurm-out-smoother_index-%j.out

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


GENOME_FILES=(genome/*.fna)
GENOME_FILE=$(echo ${GENOME_FILES[0]} | sed "s/.fna//")
INDEX_FILE=$(echo ${GENOME_FILE} | sed "s/genome/out/g")

echo "building index for ${GENOME_FILE} at resolution ${RESO}, subsamples ${SUBSAMPLE_FRACTION}, and extra_params ${EXTRA_SMOOTHER_PARAMS}"

INPUT_FOLDER=fasta
OUTPUT_FOLDER=data

TIME_FILE=out/2_index_build_times
STATUS_FILE=out/index_build_status
INDEX_SIZE_FILe=out/index_size

EXTRA_SMOOTHER_PARAMS_FILE=$(echo ${EXTRA_SMOOTHER_PARAMS} | sed "s/ /./g")

FULL_SUFFIX=${RESO}.${SUBSAMPLE_FRACTION}${EXTRA_SMOOTHER_PARAMS_FILE}

# reset file
echo -n "" > ${TIME_FILE}.${FULL_SUFFIX}.txt
echo -n "Error or in progress" > ${STATUS_FILE}.${FULL_SUFFIX}.txt
echo -n "init:" >> ${TIME_FILE}.${FULL_SUFFIX}.txt
rm -r ${INDEX_FILE}.${FULL_SUFFIX}.smoother_index

# exit on error from here on
set -e

{ time biosmoother init -d ${RESO} ${INDEX_FILE}.${FULL_SUFFIX} ${GENOME_FILE}.sizes \
                       ${GENOME_FILE}.gff 2>&1; } 2>> ${TIME_FILE}.${FULL_SUFFIX}.txt

for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
do
    SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
   echo -n "repl ${SAMPLE_NAME}:" >> ${TIME_FILE}.${FULL_SUFFIX}.txt
   { time ( zcat ${OUTPUT_FOLDER}/${SAMPLE_NAME}.${SUBSAMPLE_FRACTION}.pairs.gz | \
               biosmoother repl ${EXTRA_SMOOTHER_PARAMS} ${INDEX_FILE}.${FULL_SUFFIX} - ${SAMPLE_NAME} 2>&1 ); } 2>> \
       ${TIME_FILE}.${FULL_SUFFIX}.txt
done

du -s ${INDEX_FILE}.${FULL_SUFFIX}.smoother_index > ${INDEX_SIZE_FILe}.${FULL_SUFFIX}.txt
du -hs ${INDEX_FILE}.${FULL_SUFFIX}.smoother_index >> ${INDEX_SIZE_FILe}.${FULL_SUFFIX}.txt

echo "done"
echo -n "OK" > ${STATUS_FILE}.${FULL_SUFFIX}.txt

# biosmoother init -d 10000 HGAP3_Tb427v10.10000.even_1.smoother_index_2 ../genome/HGAP3_Tb427v10.sizes_2 ../genome/HGAP3_Tb427v10.gff_2
# cat ../data/SRR7721317.even.pairs.gz_2 | biosmoother repl HGAP3_Tb427v10.10000.even_1.smoother_index_2.smoother_index/ - R1
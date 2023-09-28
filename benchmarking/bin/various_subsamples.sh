#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=subsampling
#SBATCH --mem=250G
#SBATCH --time=12-00:00:00
#SBATCH -o out/slurm-out-subsampling-%j.out

### Pipeline according to open2C/distiller

module load ngs/bwa/0.7.16
module load ngs/samtools/1.9

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


### SET VARAIBLES ###
GENOME_FILES=(genome/*.fna)
GENOME_FILE=$(echo ${GENOME_FILES[0]} | sed "s/.fna//")

INPUT_FOLDER=fasta
OUTPUT_FOLDER=data
INT_CNT_FILE=out/interaction_count


for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
do
    READ_TWO=$(echo ${READ_FILE} | sed "s/_R1/_R2/")
    SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
    for SUBSAMPLE_FACTION in 0.2 0.4 0.6 0.8 1 #0.1 0.01 0.001 0.0001
    do
        if [ ! -f ${OUTPUT_FOLDER}/${SAMPLE_NAME}.${SUBSAMPLE_FACTION}.pairs.gz ]
        then
            bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.pairsam.gz | \
                pairtools sample ${SUBSAMPLE_FACTION} | \
                pairtools split --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}.${SUBSAMPLE_FACTION}.pairs.gz
        fi

        AVG_COV=$(bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.${SUBSAMPLE_FACTION}.pairs.gz | \
                python ../../bin/unique_interactions_per_sq_kb.py - ${GENOME_FILE}.sizes)
        echo "unique interactions per square kbp after subsampling to ${SUBSAMPLE_FACTION} = ${AVG_COV}"
        gunzip -c ${OUTPUT_FOLDER}/${SAMPLE_NAME}.${SUBSAMPLE_FACTION}.pairs.gz | \
            python ../../bin/count_unique_interactions.py - > ${INT_CNT_FILE}.${SUBSAMPLE_FACTION}.txt
    done
done

echo "done"

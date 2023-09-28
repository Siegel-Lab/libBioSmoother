#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=distiller
#SBATCH --mem=250G
#SBATCH --time=12-00:00:00
#SBATCH --partition=slim18
#SBATCH -o out/slurm-out-distiller-%j.out

### Pipeline according to open2C/distiller

module load ngs/bwa/0.7.16
module load ngs/samtools/1.9
module load ngs/sratoolkit/2.10.0

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


### SET VARAIBLES ###
GENOME_FILES=(genome/*.fna)
GENOME_FILE=$(echo ${GENOME_FILES[0]} | sed "s/.fna//")

INPUT_FOLDER=fasta
OUTPUT_FOLDER=data
INT_CNT_FILE=out/interaction_count
COVERAGE_FILE=out/unique_interactions_per_sq_kb.txt
COVERAGE_FILE_SUB=out/unique_interactions_per_sq_kb_post_subsample

NCORES=16

TARGET_COVERAGE_1=0.001
TARGET_COVERAGE_2=0.00001


### INDEX GENOME ### 
if [ ! -f ${GENOME_FILE}.fna.pac ]
then
    echo "building index"
    bwa index ${GENOME_FILE}.fna
fi

if [ ! -f ${GENOME_FILE}.sizes ]
then
    echo "counting contig sizes"
    faidx ${GENOME_FILE}.fna -i chromsizes > ${GENOME_FILE}.sizes
fi

### ALIGN READS ###
for READ_FILE in $(ls ${INPUT_FOLDER} | grep _R1)
do
    READ_TWO=$(echo ${READ_FILE} | sed "s/_R1/_R2/")
    SAMPLE_NAME=$(echo ${READ_FILE} | sed "s/_R1.fq.gz//")
    echo "working on ${SAMPLE_NAME}"
    if [ ! -f ${OUTPUT_FOLDER}/${SAMPLE_NAME}.pairsam.gz ]
    then
        echo "generating pairsam"
        TASK_TMP_DIR=$(mktemp -d -p ./ distiller.tmp.XXXXXXXXXX) 
        bwa mem -t ${NCORES} -v 3 -SP ${GENOME_FILE}.fna ${INPUT_FOLDER}/${READ_FILE} ${INPUT_FOLDER}/${READ_TWO} | \
            pairtools parse --drop-sam --min-mapq 0 --add-columns mapq,XA --walks-policy mask \
                            -c ${GENOME_FILE}.sizes | \
            pairtools select '(pair_type!="WW") and (pair_type!="XX") and not wildcard_match(pair_type, "N*")' | \
            pairtools sort --nproc ${NCORES} --tmpdir $TASK_TMP_DIR | \
            pairtools dedup -o ${OUTPUT_FOLDER}/${SAMPLE_NAME}.pairsam.gz
        rm -rf $TASK_TMP_DIR
    fi

    AVG_COV=$(bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.pairsam.gz | \
            python ../../bin/unique_interactions_per_sq_kb.py - ${GENOME_FILE}.sizes)
    echo "unique interactions per square kbp = ${AVG_COV}"
    echo "${AVG_COV}" > ${COVERAGE_FILE}
    SUBSAMPLE_FACTION_1=$(echo "scale=20 ; $TARGET_COVERAGE_1 / $AVG_COV" | bc)
    echo "fraction of pairs to use for target 1 = ${SUBSAMPLE_FACTION_1}"

    if [ ! -f ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_1.pairs.gz ]
    then
        echo "generating even_1.pairs"
        bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.pairsam.gz | \
            pairtools sample 0${SUBSAMPLE_FACTION_1} | \
            pairtools split --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_1.pairs.gz

    fi
    bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_1.pairs.gz | \
        python ../../bin/count_unique_interactions.py - > ${INT_CNT_FILE}.even_1.txt

    AVG_COV_1=$(bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_1.pairs.gz | \
            python ../../bin/unique_interactions_per_sq_kb.py - ${GENOME_FILE}.sizes)
    echo "unique interactions per square kbp after subsampling = ${AVG_COV_1} (target 1)"
    echo "target 1 is ${TARGET_COVERAGE_1}"
    echo "${AVG_COV_1}" > ${COVERAGE_FILE_SUB}_1.txt


    SUBSAMPLE_FACTION_2=$(echo "scale=20 ; $TARGET_COVERAGE_2 / $AVG_COV" | bc)
    echo "fraction of pairs to use for target 1 = ${SUBSAMPLE_FACTION_2}"

    if [ ! -f ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_2.pairs.gz ]
    then
        echo "generating even_2.pairs"
        bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.pairsam.gz | \
            pairtools sample 0${SUBSAMPLE_FACTION_2} | \
            pairtools split --output-pairs ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_2.pairs.gz

    fi
    bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_2.pairs.gz | \
        python ../../bin/count_unique_interactions.py - > ${INT_CNT_FILE}.even_2.txt

    AVG_COV_2=$(bgzip -cd -@ 3 ${OUTPUT_FOLDER}/${SAMPLE_NAME}.even_2.pairs.gz | \
            python ../../bin/unique_interactions_per_sq_kb.py - ${GENOME_FILE}.sizes)
    echo "unique interactions per square kbp after subsampling = ${AVG_COV_2} (target 2)"
    echo "target 1 is ${TARGET_COVERAGE_2}"
    echo "${AVG_COV_2}" > ${COVERAGE_FILE_SUB}_2.txt
done

echo "done"

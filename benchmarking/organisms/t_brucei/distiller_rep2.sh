#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --mail-user=markus.rainer.schmidt@gmail.com
#SBATCH --mail-type=END
#SBATCH --job-name=distiller
#SBATCH --mem=250G
#SBATCH --time=12-00:00:00
#SBATCH -o out/slurm-out-distiller-%j.out

module load ngs/bwa/0.7.16
module load ngs/samtools/1.9
module load ngs/sratoolkit/2.10.0

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother

NCORES=16

GENOME_FILES=(genome/*.fna)
GENOME_FILE=$(echo ${GENOME_FILES[0]} | sed "s/.fna//")
OUTPUT_FOLDER=data

TASK_TMP_DIR=$(mktemp -d -p ./ distiller.tmp.XXXXXXXXXX) 

bwa mem -t ${NCORES} -v 3 -SP ${GENOME_FILE}.fna fasta/rep2/SRR7721307_R1.fq.gz fasta/rep2/SRR7721307_R2.fq.gz | \
    pairtools parse --drop-sam --min-mapq 0 --add-columns mapq,XA --walks-policy mask \
                    -c ${GENOME_FILE}.sizes | \
    pairtools select '(pair_type!="WW") and (pair_type!="XX") and not wildcard_match(pair_type, "N*")' | \
    pairtools sort --nproc ${NCORES} --tmpdir $TASK_TMP_DIR | \
    pairtools dedup | \
    pairtools split --output-pairs ${OUTPUT_FOLDER}/SRR7721307.1.pairs.gz

rm -rf $TASK_TMP_DIR
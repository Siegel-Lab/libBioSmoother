#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=smoother_index
#SBATCH --mem=250G
#SBATCH --partition=slim16
#SBATCH --time=12-00:00:00
#SBATCH -o out/slurm-out-smoother_index-%j.out

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


GENOME_FILES=(genome/*.fna)
GENOME_FILE=$(echo ${GENOME_FILES[0]} | sed "s/.fna//")
INDEX_FILE=out/comparison_index
RESO=10000


# reset file
rm -r ${INDEX_FILE}.smoother_index

biosmoother init ${INDEX_FILE} ${GENOME_FILE}.sizes ${GENOME_FILE}.gff -d ${RESO}

echo "a"
zcat data/merged.1.pairs.gz | biosmoother repl ${INDEX_FILE} - WT
echo "b"
zcat data/SRR7721307.1.pairs.gz | biosmoother repl ${INDEX_FILE} - H3.V-/-
echo "c"


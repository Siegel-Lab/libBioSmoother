#!/bin/bash

#source ~/.miniconda3/etc/profile.d/conda.sh
#conda activate smoother

#module load ngs/sratoolkit/2.10.0

# tryps
# Paper https://www.nature.com/articles/s41586-018-0619-8#data-availability
# Ncbi https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA393359&o=assay_type_s%3Aa

#mkdir fasta
cd fasta


#fasterq-dump SRR7721318   # WT
#fasterq-dump SRR7721304   # TetR, T7RNAP
#fasterq-dump SRR7721305   # TetR, T7RNAP
#fasterq-dump SRR7721306   # TetR, T7RNAP

for R in *.fastq; do
    echo $R
    gzip $R &
done

wait

# merged is SRR7721317 and SRR7721318 together (both WT)
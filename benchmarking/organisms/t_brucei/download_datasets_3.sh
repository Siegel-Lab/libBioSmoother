#!/bin/bash

#source ~/.miniconda3/etc/profile.d/conda.sh
#conda activate smoother

module load ngs/sratoolkit/2.10.0

# tryps
# Paper https://www.nature.com/articles/s41586-018-0619-8#data-availability
# Ncbi https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA393359&o=assay_type_s%3Aa

#mkdir fasta
cd fasta


fasterq-dump SRR7721307  # TetR, T7RNAP, H3.V-/-

# for R in *.fastq; do
#     echo $R
#     gzip $R &
# done

# wait

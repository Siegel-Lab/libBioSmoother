#!/bin/bash

# tryps
# Paper https://www.nature.com/articles/s41586-018-0619-8#data-availability
# Ncbi https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA393359&o=assay_type_s%3Aa
module load ngs/sratoolkit/2.10.0
# fasterq-dump SRR7988175 # WT rna-seq
# mv SRR7988175_1.fastq fasta/disabled/SRR7988175_R1.fq
# mv SRR7988175_2.fastq fasta/disabled/SRR7988175_R2.fq
# srun gzip fasta/disabled/SRR7988175_R1.fq
srun gzip fasta/disabled/SRR7988175_R2.fq


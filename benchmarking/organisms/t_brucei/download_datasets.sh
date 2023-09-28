#!/bin/bash

# tryps
# Paper https://www.nature.com/articles/s41586-018-0619-8#data-availability
# Ncbi https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA393359&o=assay_type_s%3Aa

mkdir fasta
cd fasta
# 26 SRR7721317 PRJNA486752 SAMN09863016 Hi-C 14.05 G 4.66 Gb SRX4578221 GSM3346686 other GENOMIC Trypanosoma brucei brucei 2018-09-25 2018-08-2011:16:00Z GSM3346686 whole liquid cell culture SRP158382 Bloodstream form Lister 427 MITat 1.2 brucei wild-type

fasterq-dump SRR7721317
mv SRR7721317_1.fastq SRR7721317_1.fq
mv SRR7721317_2.fastq SRR7721317_2.fq
srun gzip *.fq
cd ..

mkdir genome


#fasterq-dump SRR7721318    WT
#fasterq-dump SRR7721304    TetR, T7RNAP
#fasterq-dump SRR7721305    TetR, T7RNAP
#fasterq-dump SRR7721306    TetR, T7RNAP
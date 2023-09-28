#!/bin/bash


source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother

module load ngs/sratoolkit/2.10.0

# e. coli
# Paper https://www.sciencedirect.com/science/article/pii/S0092867417315076
# ncbi https://www.ncbi.nlm.nih.gov/sra/SRX3451241[accn]

cd fasta
# accession: GSM2870439
fasterq-dump SRR6354577
srun gzip *.fq
cd ..


cd genome

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
#srun gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
#srun faidx GCF_000005845.2_ASM584v2_genomic.fna -i chromsizes > GCF_000005845.2_ASM584v2_genomic.sizes

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
#srun gunzip GCF_000005845.2_ASM584v2_genomic.gff.gz
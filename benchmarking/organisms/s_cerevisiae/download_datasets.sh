#!/bin/bash

# YEAST
# Paper https://www.nature.com/articles/s41467-021-26629-6
# Ncbi https://www.ncbi.nlm.nih.gov/sra?term=SRP284851

mkdir fasta
cd fasta
# SRX10818035: GSM5287083: Sth1-K501R_G1+IAA_rep2; Saccharomyces cerevisiae; Hi-C
fasterq-dump SRR14469367
gzip *.fq
cd ..

mkdir genome
cd genome
# genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
gunzip GCF_000146045.2_R64_genomic.fna.gz
# genome sizes
faidx GCF_000146045.2_R64_genomic.fna -i chromsizes > GCF_000146045.2_R64_genomic.sizes

# annotation
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz
gunzip GCF_000146045.2_R64_genomic.gff.gz 
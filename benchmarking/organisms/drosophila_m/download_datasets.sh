#!/bin/bash


source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother

module load ngs/sratoolkit/2.10.0

# drosophila
# Paper https://www.nature.com/articles/s41467-017-02526-9
# ncbi https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE101317&o=acc_s%3Aa

mkdir fasta
cd fasta
# accession: GSE101317
fasterq-dump SRR5820092
srun gzip *.fq
cd ..

cd ..

#mkdir genome
#cd genome

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
#srun gunzip GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz
#srun faidx GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna -i chromsizes > GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.sizes

#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz
#srun gunzip GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff.gz
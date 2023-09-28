#!/bin/bash


source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother

module load ngs/sratoolkit/2.10.0

# human
# Paper https://www.cell.com/fulltext/S0092-8674(14)01497-4#secsectitle0140
# ncbi https://www.ncbi.nlm.nih.gov/Traces/study/?acc=GSE63525&o=acc_s%3Aa

cd fasta
# accession: GSE63525
#fasterq-dump SRR1658572
srun gzip *.fq
cd ..

exit

cd genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

srun gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
srun faidx GCF_000001405.40_GRCh38.p14_genomic.fna -i chromsizes > GCF_000001405.40_GRCh38.p14_genomic.size

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz
srun gunzip GCF_000001405.40_GRCh38.p14_genomic.gff.gz
#!/bin/bash
#SBATCH --cpus-per-task=18
#SBATCH --job-name=smoother_index
#SBATCH --mem=250G
#SBATCH --partition=slim16
#SBATCH --time=12-00:00:00
#SBATCH -o slurm-out-smoother_index-%j.out

source ~/.miniconda3/etc/profile.d/conda.sh
conda activate smoother


## copy files 

#ANNA_FOLDER=/work/project/ladsie_012/2023_Smoother/Bonetti2020_RADICL
# cp ${ANNA_FOLDER}/Bonetti2020_RADICL.PRE2 Bonetti2020_RADICL.PRE2
# cp ${ANNA_FOLDER}/GCF_000001635.27_GRCm39_genomic.gff GCF_000001635.27_GRCm39_genomic.gff
# cp ${ANNA_FOLDER}/GCF_000001635.27_GRCm39_genomic.fna GCF_000001635.27_GRCm39_genomic.fna
# faidx GCF_000001635.27_GRCm39_genomic.fna -i chromsizes > GCF_000001635.27_GRCm39_genomic.sizes

rm -r index.50000.smoother_index
rm -r index.10000.smoother_index

# create index
biosmoother init -d 50000 index.50000 GCF_000001635.27_GRCm39_genomic.sizes GCF_000001635.27_GRCm39_genomic.gff
biosmoother repl index.50000.smoother_index Bonetti2020_RADICL.PRE2 Bonetti2020 -C readid chr1 pos1 chr2 pos2 mapq1 mapq2 xa1 xa2

biosmoother init -d 10000 index.10000 GCF_000001635.27_GRCm39_genomic.sizes GCF_000001635.27_GRCm39_genomic.gff
biosmoother repl index.10000.smoother_index Bonetti2020_RADICL.PRE2 Bonetti2020 -C readid chr1 pos1 chr2 pos2 mapq1 mapq2 xa1 xa2
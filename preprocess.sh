#!/bin/bash
#SBATCH --mem 100G -J preprocess_smoother --time=240:00:00 -o slurm_preprocess_heatmap-%j.out --mail-user=markus.rainer.schmidt@gmail.com --mail-type END

if true;
then
    # source activate $(pwd)/conda_env/smoother
    source activate libContactMapping

    #./bin/conf_version.sh
    #cat VERSION

    BEDS="/work/project/ladsie_012/ABS.2.2/2021-10-26_NS502-NS521_ABS_CR_RADICL_inputMicroC/bed_files/minus_N"
    BED_SUF="RNA.PRE_K1K2.PRE_K3.PRE_R_D.PRE_R_D_K1K2.PRE_R_D_PRE2.bedsorted.PRE2"

    BAMS="/work/project/ladsie_012/ABS.2.2/20210608_Inputs"
    BAM_SUF="R1.sorted.bam"

    INDEX_PREFIX="../smoother_out/radicl-2"

    rm -r ${INDEX_PREFIX}.smoother_index


    echo "working on index ${INDEX_PREFIX}"


    #python3 python/main.py indexer init "${INDEX_PREFIX}" "../smoother_in/Lister427.sizes" -d 1000
    python3 python/main.py indexer init "${INDEX_PREFIX}" "../smoother_in/Lister427_no_uni.sizes" -d 10


    python3 python/main.py indexer anno "${INDEX_PREFIX}" "../smoother_in/HGAP3_Tb427v10_merged_2021_06_21.gff3"

    #gdb python3 -ex "run python/main.py indexer repl \"${INDEX_PREFIX}\" \"${BEDS}/NS503_P10_Total_2.${BED_SUF}\" \"P10_Total_Rep2\" -g a"
    python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS503_P10_Total_2.${BED_SUF}" "P10_Total_Rep2" -g a

    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS617_P10_Total_1.${BED_SUF}" "P10_Total_Rep1" -g a
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS504_P10_Total_3.${BED_SUF}" "P10_Total_Rep3" -g a
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS505_N50_Total_1.${BED_SUF}" "N50_Total_Rep1" -g a
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS506_N50_Total_2.${BED_SUF}" "N50_Total_Rep2" -g a
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS507_N50_Total_3.${BED_SUF}" "N50_Total_Rep3" -g a
    #
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS508_P10_NPM_1.${BED_SUF}" "P10_NPM_Rep1" -g b
    #
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS509_P10_NPM_2.${BED_SUF}" "P10_NPM_Rep2" -g b
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS510_P10_NPM_3.${BED_SUF}" "P10_NPM_Rep3" -g b
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS511_N50_NPM_1.${BED_SUF}" "N50_NPM_Rep1" -g b
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS512_N50_NPM_2.${BED_SUF}" "N50_NPM_Rep2" -g b
    #python3  python/main.py indexer repl "${INDEX_PREFIX}" "${BEDS}/NS513_N50_NPM_3.${BED_SUF}" "N50_NPM_Rep3" -g b

else
    rm -r ../smoother_out/hic2.smoother_index

    python3 python/main.py indexer init ../smoother_out/hic2 ../smoother/Lister427.sizes -d 1000 --test

    python3 python/main.py indexer anno ../smoother_out/hic2 ../smoother/HGAP3_Tb427v10_merged_2021_06_21.gff3

    #python3 python/main.py indexer repl ../smoother_out/radicl ../smoother_in/anna.sort.test.PRE2 P10_Total
    #gdb python3 -ex "run python/main.py indexer repl ../smoother_out/hic2 ../smoother_in/claudia.pre1 P10_R1"

    python3 python/main.py indexer repl -q -m ../smoother_out/hic2 ../smoother_in/claudia.pre1 P10_R1


    #gdb python3 -ex "run python/main.py indexer track ../smoother_out/hic ../smoother_in/coverage.tsv.sorted rna_seq"
    #python3 python/main.py indexer track ../smoother_out/hic ../smoother_in/coverage.tsv.sorted rna_seq
fi
#!/bin/bash

# while [[ $(squeue | grep $(whoami)) ]]
# do
#     echo "jobs running, waiting 60 seconds"
#     sleep 60
# done
# echo "no jobs running, starting benchmarking"

for RESO in 50000 10000 5000
do
    for SUBSAMPLE_FRACTION in "even_1"
    do
        for FOLDER in drosophila_m e_coli t_brucei s_cerevisiae
        do
            # try 4 configurations:
            #
            # ID:               1   2   3   4
            # mapping_q         +   +   -   -
            # multimappers      +   -   +   -
            # category          +   -   -   -
            # 1: ""
            # 2: "-c -m"
            # 3: "-q -c"
            # 4: "-q -c -m"

            cd organisms/${FOLDER}
            for EXTRA_SMOOTHER_PARAMS in "" #" -c -m" " -q -c" " -q -c -m"
            do
                echo "benchmarking index for ${FOLDER} at resolution ${RESO} with extra parameters ${EXTRA_SMOOTHER_PARAMS}"
                sbatch --export ALL,SUBSAMPLE_FRACTION=${SUBSAMPLE_FRACTION},RESO=${RESO},EXTRA_SMOOTHER_PARAMS="${EXTRA_SMOOTHER_PARAMS}" ../../bin/benchmark_index.sh 
            done
            cd ../..
        done
        
        cd organisms/t_brucei
        for Q in " -q" ""
        do
            for C in " -c" ""
            do
                for M in " -m" ""
                do
                    EXTRA_SMOOTHER_PARAMS="${Q}${C}${M}"
                    # two of these have already been triggered above
                    if [ "$EXTRA_SMOOTHER_PARAMS" != "" ] # && \
                    # [ "$EXTRA_SMOOTHER_PARAMS" != " -c -m" ] && \
                    # [ "$EXTRA_SMOOTHER_PARAMS" != " -q -c" ] && \
                    # [ "$EXTRA_SMOOTHER_PARAMS" != " -q -c -m" ]
                    then
                        echo "benchmarking index at resolution ${RESO} and extra params ${EXTRA_SMOOTHER_PARAMS}"
                        sbatch --export ALL,SUBSAMPLE_FRACTION=${SUBSAMPLE_FRACTION},RESO=${RESO},EXTRA_SMOOTHER_PARAMS="${EXTRA_SMOOTHER_PARAMS}" ../../bin/benchmark_index.sh 
                    fi
                done
            done
        done
        cd ../..
    done


    cd organisms/t_brucei
    for SUBSAMPLE_FRACTION in 0.2 0.4 0.6 0.8 1
    do
        for EXTRA_SMOOTHER_PARAMS in ""
        do
            echo "benchmarking index at resolution ${RESO} and subsamples ${SUBSAMPLE_FRACTION} at extra params ${EXTRA_SMOOTHER_PARAMS}"
            sbatch --export ALL,SUBSAMPLE_FRACTION=${SUBSAMPLE_FRACTION},RESO=${RESO},EXTRA_SMOOTHER_PARAMS="${EXTRA_SMOOTHER_PARAMS}" ../../bin/benchmark_index.sh
        done
    done
    cd ../.. 
done

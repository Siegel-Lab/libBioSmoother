#!/bin/bash


for FOLDER in $(ls organisms)
do
    cd organisms/${FOLDER}
    echo "running distiller on ${FOLDER}"
    sbatch ../../bin/distiller.sh
    cd ../..
done
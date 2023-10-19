#!/bin/bash

# run this with source ./publish_mac_m1.sh

set -e
if [ ! -d "libSps" ]; then
    git clone https://github.com/Siegel-Lab/libSps.git
fi

for PY_VERSION in "3.9" "3.10" "3.11"
do 

    conda create -n tmp_publish_smoother python=$PY_VERSION -y
    conda activate tmp_publish_smoother

    python -m pip install --upgrade pip
    python -m pip install --upgrade wheel pybind11 cmake

    cd libSps
    pip install -e .
    cd ..

    python3 setup.py bdist_wheel

    conda deactivate
    conda env remove -n tmp_publish_smoother -y

done

conda create -n tmp_publish_smoother python=3.9 -y
conda activate tmp_publish_smoother

python -m pip install --upgrade wheel twine

twine upload dist/*

conda deactivate
conda env remove -n tmp_publish_smoother -y
    
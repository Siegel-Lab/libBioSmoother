#!/bin/bash

# run this with source ./test_publish_mac_m1.sh

set -e

for PY_VERSION in "3.8" "3.9" "3.10" "3.11"
do 

    conda create -n tmp_publish_smoother python=$PY_VERSION -y
    conda activate tmp_publish_smoother

    python -m pip install --upgrade pip
    python -m pip install --upgrade wheel pybind11 cmake

    git clone https://github.com/Siegel-Lab/libSps.git
    cd libSps
    pip install -e .
    cd ..

    python3 setup.py bdist_wheel

    conda deactivate
    conda env remove -n tmp_publish_smoother -y
done

conda create -n tmp_publish_smoother python=3.8 -y
conda activate tmp_publish_smoother

python -m pip install --upgrade wheel twine

twine upload -r testpypi dist/*

conda deactivate
conda env remove -n tmp_publish_smoother -y

for PY_VERSION in "3.8" "3.9" "3.10" "3.11"
do 

    conda create -n tmp_test_smoother python=$PY_VERSION -y
    conda activate tmp_test_smoother

    python3 -m pip install scipy statsmodles scikit-learn drawSvg==1.9.0

    python3 -m pip install --index-url https://test.pypi.org/simple/ libbiosmoother --no-binary :all:

    libbiosmoother --help

    conda deactivate
    conda env remove -n tmp_publish_smoother -y

done

echo "All tests passed!"
#!/bin/bash

# build conda package
conda build . -c conda-forge

conda install --use-local libsmoother
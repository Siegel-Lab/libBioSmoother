#!/bin/bash

# Compile C++ code

# create build directory
mkdir -p build
cd build

# compile code
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make VERBOSE=1 -j 8

# move libary into python package
cp libSmoother.*.so ../libSmoother/
cd ..

# install python package
python setup.py install
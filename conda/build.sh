#!/bin/bash

# Compile C++ code

# create build directory
mkdir -p build
cd build

# compile code
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make VERBOSE=1 -j 8

# move libary into python package
mkdir -p $PREFIX/lib
cp libsmoothercpp*.so $PREFIX/lib
ln -s $PREFIX/lib/libsmoothercpp*.so $PREFIX/lib/libsmoothercpp.so

# install python package
cd ..
$PYTHON conda/setup.py install
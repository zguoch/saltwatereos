#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
set -x

rootpath=${PWD}
export GSL_ROOT_DIR=${rootpath}/gsl-2.6/gsl
# 1. build gsl
build_gsl=$rootpath/gsl-2.6/src/build 
mkdir $build_gsl
cd $build_gsl 
../configure --prefix=${GSL_ROOT_DIR}
make -j 8 #using 8 threads
make install

# 2. build swEOS
build_swEOS=$rootpath/build 
mkdir $build_swEOS
cd $build_swEOS
cmake ..
make install
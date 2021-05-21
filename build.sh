#!/bin/env bash

##################
# Set the compiler
##################

# For GNU
export CC=gcc
export FC=gfortran

# For Intel
#export CC=icc
#export FC=ifort

####################
# Set the build type
####################

export BUILD_TYPE=debug
#export BUILD_TYPE=release

#######################
# Clean build directory
#######################
rm -rf [B-Zc-z]*

###########
# Run cmake
###########
cmake .. -DCMAKE_BUILD_TYPE=${BUILD_TYPE}

#######
# Build
#######
make -j4 VERBOSE=1

######
# Test
######
ctest --output-on-failure

#!/bin/bash
PYDOF_HOME=$PWD
PYDOF_SRC=$PYDOF_HOME/src
PYDOF_EXAMPLES=$PYDOF_HOME/examples/
PYDOF_FSRC=$PYDOF_SRC/fortran/
PYDOF_CSRC=$PYDOF_SRC/c++/
FCOMP='gfortran -fopenmp -O3 -mcmodel=large -fcheck=all'
CXX='g++'
CFLAGS="-O3"
INCLUDES="-I$PYDOF_CSRC"

echo ''
echo ' -------------------'
echo ' PYDOF'
echo ' Author: Z. Nikolaou'
echo ' -------------------'
echo ''
echo 'PYDOF home dir:'
echo ' '$PYDOF_HOME
echo 'PYDOF source:'
echo ' '$PYDOF_SRC
echo 'PYDOF Fortran source:'
echo ' '$PYDOF_FSRC
echo 'PYDOF C++ source:'
echo ' '$PYDOF_CSRC
echo 'PYDOF examples:'
echo ' '$PYDOF_EXAMPLES
echo 'Fortran compiler:'
echo ' '$FCOMP
echo 'C++ compiler:'
echo ' '$CXX
echo ''
echo 'Environment variables set! '
echo ''

export PYDOF
export PYDOF_HOME
export PYDOF_SRC
export PYDOF_FSRC
export PYDOF_CSRC
export PYDOF_EXAMPLES
export FCOMP
export CXX
export CFLAGS
export INCLUDES

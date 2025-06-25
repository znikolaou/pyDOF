#!/bin/bash
PYDOF_HOME=$PWD
PYDOF_SRC=$PYDOF_HOME/src
PYDOF_TUTORIALS=$PYDOF_HOME/tutorials/
PYDOF_FSRC=$PYDOF_SRC/fortran/
PYDOF_CSRC=$PYDOF_SRC/c++/
FCOMP='gfortran -fopenmp -O3 -mcmodel=large -fcheck=all'
CXX='g++'
CFLAGS="-O3"
INCLUDES="-I$PYDOF_CSRC"

echo 'PYDOF home dir:'
echo $PYDOF_HOME
echo 'PYDOF source:'
echo $PYDOF_SRC
echo 'PYDOF Fortran source:'
echo $PYDOF_FSRC
echo 'PYDOF C++ source:'
echo $PYDOF_CSRC
echo 'PYDOF tutorials:'
echo $PYDOF_TUTORIALS
echo 'Fortran compiler:'
echo $FCOMP
echo 'C++ compiler:'
echo $CXX

export PYDOF
export PYDOF_HOME
export PYDOF_SRC
export PYDOF_FSRC
export PYDOF_CSRC
export PYDOF_TUTORIALS
export FCOMP
export CXX
export CFLAGS
export INCLUDES

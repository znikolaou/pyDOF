#!/bin/bash
PYDOF_HOME=$PWD
PYDOF_SRC=$PYDOF_HOME/src
PYDOF_TUTORIALS=$PYDOF_HOME/tutorials/
PYDOF_FSRC=$PYDOF_SRC/fortran/
FCOMP='gfortran -fopenmp -O3 -mcmodel=large -fcheck=all'

echo 'PYDOF home dir:'
echo $PYDOF
echo 'PYDOF source:'
echo $PYDOF_SRC
echo 'PYDOF Fortran source:'
echo $PYDOF_FSRC
echo 'PYDOF tutorials:'
echo $PYDOF_TUTORIALS
echo 'Fortran compiler'
echo $FCOMP

export PYDOF
export PYDOF_HOME
export PYDOF_SRC
export PYDOF_FSRC
export PYDOF_TUTORIALS
export FCOMP

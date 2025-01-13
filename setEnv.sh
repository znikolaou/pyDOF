#!/bin/bash
PYDOF=$PWD
PYDOF_SRC=$PYDOF/src
PYDOF_FORTRAN_SRC=$PYDOF_SRC/fortran/
PYDOF_TUTORIALS=$PYDOF/tutorials/
#
PYDOF_FCOMPILER=gfortran 
PYDOF_FOPTIONS='-O3 -mcmodel=large'

echo 'PYDOF:'
echo $PYDOF
echo 'PYDOF source:'
echo $PYDOF_SRC
echo 'PYDOF Fortran source:'
echo $PYDOF_FORTRAN_SRC
echo 'PYDOF tutorials:'
echo $PYDOF_TUTORIALS
echo 'Fortran compiler'
echo $PYDOF_FCOMPILER
echo 'Fortran options:'
echo $PYDOF_FOPTIONS

export PYDOF
export PYDOF_SRC
export PYDOF_FORTRAN_SRC
export PYDOF_TUTORIALS
export PYDOF_FCOMPILER
export PYDOP_FOPTIONS

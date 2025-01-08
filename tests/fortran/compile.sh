#!/bin/bash

FC=gfortran
OPTS=-O3

# $FC $OPTS testReadCoefs.f $PYDOF_FORTRAN_SRC/pydofParser.f -o goTestReadCoeffs

$FC $OPTS testChecks.f $PYDOF_FORTRAN_SRC/pydof_filter.f -o goTestChecks

$FC $OPTS testReadCoefs.f $PYDOF_FORTRAN_SRC/pydof_filter.f $PYDOF_FORTRAN_SRC/pydof_parser.f -o goTestReadCoeffs



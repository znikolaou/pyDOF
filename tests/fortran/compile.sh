#!/bin/bash

FC=gfortran
OPTS=-O3

# $FC $OPTS testReadCoefs.f $FORTRAN_SRC/pydofParser.f -o goTestReadCoeffs

$FC $OPTS testChecks.f $FORTRAN_SRC/pydof_filter.f -o goTestChecks

$FC $OPTS testReadCoefs.f $FORTRAN_SRC/pydof_filter.f $FORTRAN_SRC/pydof_parser.f -o goTestReadCoeffs



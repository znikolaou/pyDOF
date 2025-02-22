#!/bin/bash

rm go* log* data*

#$PYDOF_FCOMPILER $PYDOF_FOPTIONS -o goTestHeaderFile testHeaderFile.f 
#$PYDOF_FCOMPILER -o goTestBcs testBcs.f $PYDOF_FORTRAN_SRC/pydof_filter.f -I/$PWD
#$PYDOF_FCOMPILER $PYDOF_FOPTIONS-o goTestWrkSize testWrkSize.f $PYDOF_FORTRAN_SRC/pydof_filter.f -I/$PWD
gfortran -fopenmp -o goTestFilter1d testFilter1d.f $PYDOF_FORTRAN_SRC/pydof_filter.f -I/$PWD


#$PYDOF_FCOMPILER $PYDOF_FOPTIONS -o goTestFilter1d testFilter1d.f $PYDOF_FORTRAN_SRC/pydof_filter.f -I/$PWD





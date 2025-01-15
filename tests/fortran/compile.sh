#!/bin/bash

$PYDOF_FCOMPILER $PYDOF_FOPTIONS -o goTestChecks testChecks.f $PYDOF_FORTRAN_SRC/pydof_filter.f -I/$PWD
$PYDOF_FCOMPILER $PYDOF_FOPTIONS -o goTestHeaderFile testHeaderFile.f 
$PYDOF_FCOMPILER -o goTestBcs testBcs.f $PYDOF_FORTRAN_SRC/pydof_filter.f -I/$PWD




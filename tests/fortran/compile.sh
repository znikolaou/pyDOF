#!/bin/bash

$PYDOF_FCOMPILER $PYDOF_FOPTIONS testChecks.f $PYDOF_FORTRAN_SRC/pydof_filter.f -o goTestChecks

$PYDOF_FCOMPILER $PYDOF_FOPTIONS testHeaderFile.f -o goTestHeaderFile



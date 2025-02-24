#!/bin/bash

rm go* log* data*

$FCOMP testHeaderFile.f -o goTestHeaderFile
$FCOMP testBcs.f $PYDOF_FSRC/pydof_filter.f -I/$PWD -o goTestBcs 
$FCOMP testWrkSize.f $PYDOF_FSRC/pydof_filter.f -I/$PWD -o goTestWrkSize 
$FCOMP testFilter1d.f $PYDOF_FSRC/pydof_filter.f -I/$PWD -o goTestFilter1d 





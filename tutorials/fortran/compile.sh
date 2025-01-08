#!/bin/bash
FILTERS_LIB=./filters/pydofLib.a
gfortran main1D.f $FILTERS_LIB -o go1D

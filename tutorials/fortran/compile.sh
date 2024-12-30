#!/bin/bash
FILTERS_LIB_DIR=./filters/
gfortran main.f filters.f -I $FILTERS_LIB_DIR -o run

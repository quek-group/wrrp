#!/bin/bash

rm *mod 
gfortran -o struct_sc.x params.f90   mergesort.f90 supercell.f90 struct_sc.f90

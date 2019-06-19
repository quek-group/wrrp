#!/bin/bash
# An example script to compile the COMPLEX version of the code
# We compile with LAPACK, FFTW2, and OpenMP
# Note that you will need to modify this script to suit your machine's environment

rm *mod
ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90  -lfftw2xf_double_intel $LIBLAPACK_MT
#ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90 -I/opt/intel/mkl/include/fftw/ -L/opt/intel/mkl/lib/intel64/ -lfftw2xf_double_intel -mkl

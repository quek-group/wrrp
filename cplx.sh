#!/bin/bash
# Compile complex version of the code

 
rm *mod
#ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90 -I/opt/intel/mkl/include/fftw/ -L/opt/intel/mkl/lib/intel64/ -lfftw2xf_double_intel -mkl
ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90  -lfftw2xf_double_intel $LIBLAPACK_MT
#ifort -O3 -fpp -DCPLX -xHost -qopt-matmul  -traceback -check all -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90  -lfftw2xf_double_intel $LIBLAPACK_MT

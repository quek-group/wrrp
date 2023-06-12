#!/bin/bash
# Compile complex version of the code

rm *mod
# wrrp.x
#ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90 -I/opt/intel/mkl/include/fftw/ -L/opt/intel/mkl/lib/intel64/ -lfftw2xf_double_intel -mkl

# wrrp.sc.x
#ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wrrp.cplx.sc.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90  wrrp.f90  -lfftw2xf_double_intel $LIBLAPACK_MT

# wfcint.x
#ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o wfcint.cplx.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90 interp.f90 wfcint.f90  -I/opt/intel/mkl/include/fftw/ -L/opt/intel/mkl/lib/intel64/ -lfftw2xf_double_intel -mkl

# autocorr.x
ifort -O3 -fpp -DCPLX -xHost -qopt-matmul -qopenmp -parallel -traceback -o autocorr.cplx.x params.f90 v.f90 mergesort.f90 supercell.f90 fftw.f90 input.f90 interp.f90 autocorr.f90  -I/opt/intel/mkl/include/fftw/ -L/opt/intel/mkl/lib/intel64/ -lfftw2xf_double_intel -mkl

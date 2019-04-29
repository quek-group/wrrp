#!/bin/bash
# Compile complex version of the code

ifort -O3 -fpp -qopt-matmul -qopenmp -parallel -o wrrp.real.x wrrp.f90 input.f90 params.f90 fftw.f90 v.f90 -I/opt/intel/mkl/include/fftw/ -L/opt/intel/mkl/lib/intel64/ -lfftw2xf_double_intel -mkl

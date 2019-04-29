#!/bin/bash
# Compile complex version of the code

ifort -O3 -fpp -DCPLX -o read_datarrp.cplx.1lgr.x xsf.f90 read_datarrp.f90 

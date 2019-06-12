#!/bin/bash

WDIR=$(pwd)
PREFIX='si'

# link SCF charge density to WFN directories
mkdir -p ./2-wfn/$PREFIX.save
ln -s ./1-scf/$PREFIX.save/charge-density.hdf5   ./2-wfn/$PREFIX.save
ln -s ./1-scf/$PREFIX.save/data-file-schema.xml  ./2-wfn/$PREFIX.save

mkdir -p ./3-wfnq/$PREFIX.save
ln -s ./1-scf/$PREFIX.save/charge-density.hdf5   ./3-wfnq/$PREFIX.save
ln -s ./1-scf/$PREFIX.save/data-file-schema.xml  ./3-wfnq/$PREFIX.save

# link WFN files to EPSILON / SIGMA directories
ln -s ./2-wfn/wfn.cplx   ./4-epsilon/WFN
ln -s ./2-wfn/wfn.cplx   ./5-epsilon_qmap/WFN

ln -s ./3-wfnq/wfn.cplx  ./4-epsilon/WFNq
ln -s ./3-wfnq/wfn.cplx  ./5-epsilon_qmap/WFNq

ln -s ./2-wfn/wfn.cplx   ./6-sigma_qmap/WFN_inner

# link ancillary files to SIGMA directory
ln -s ./2-wfn/vxc.dat      ./6-sigma_qmap/vxc.dat
ln -s ./2-wfn/rho.cplx     ./6-sigma_qmap/RHO
ln -s ./4-epsilon/eps0mat  ./6-sigma_qmap/eps0mat
ln -s ./4-epsilon/epsmat   ./6-sigma_qmap/epsmat

# link files for WRRP
ln -s ./4-epsilon/eps0mat         ./7-wrrp/eps0mat
ln -s ./4-epsilon/epsmat          ./7-wrrp/epsmat
ln -s ./5-epsilon_qmap/isortg     ./7-wrrp/isortg
ln -s ./6-sigma_qmap/gmapdata     ./7-wrrp/gmapdata
ln -s ./6-sigma_qmap/wrrp.wcoul0  ./7-wrrp/wrrp.wcoul

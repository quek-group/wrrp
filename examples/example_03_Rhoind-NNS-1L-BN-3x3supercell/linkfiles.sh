#!/bin/bash

WDIR=$(pwd)
PREFIX='bn'

# link SCF charge density to WFN directories
mkdir -p ./2-wfn/$PREFIX.save
cd $WDIR/2-wfn/$PREFIX.save
ln -s ../../1-scf/$PREFIX.save/charge-density.hdf5   .
ln -s ../../1-scf/$PREFIX.save/data-file-schema.xml  .

cd $WDIR
mkdir -p ./3-wfn_nns/$PREFIX.save
cd $WDIR/3-wfn_nns/$PREFIX.save
ln -s ../../1-scf/$PREFIX.save/charge-density.hdf5   .
ln -s ../../1-scf/$PREFIX.save/data-file-schema.xml  .

# link WFN files to EPSILON / SIGMA directories
cd $WDIR/4-epsilon
ln -s ../2-wfn/wfn.cplx  ./WFN
ln -s ../3-wfn_nns/wfn.cplx ./WFNq

cd $WDIR/5-epsilon_qmap
ln -s ../2-wfn/wfn.cplx  ./WFN
ln -s ../3-wfn_nns/wfn.cplx ./WFNq

cd $WDIR/6-sigma_qmap
ln -s ../2-wfn/wfn.cplx  ./WFN_inner

# link ancillary files to SIGMA directory
ln -s ../2-wfn/vxc.dat      ./vxc.dat
ln -s ../2-wfn/rho.cplx     ./RHO
ln -s ../4-epsilon/eps0mat  ./eps0mat
ln -s ../4-epsilon/epsmat   ./epsmat

# link files for WRRP
cd $WDIR/7-wrrp
#ln -s ../2-wfn/subweights.dat  # link subweights.dat file from BGW subsample_plan.x
ln -s ../4-epsilon/eps0mat         .
ln -s ../4-epsilon/epsmat          .
ln -s ../5-epsilon_qmap/isortg     .
ln -s ../6-sigma_qmap/gmapdata     .

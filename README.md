# wrrp.x

## Description
wrrp.x is a Fortran code that performs six-dimensional (6D) inverse Fourier transforms on inverse dielectric matrices generated using [BerkeleyGW](https://berkeleygw.org/).

The BerkeleyGW/1.2.0 code is licensed under licensed under a free, open source, and permissive 3-clause modified BSD license. BerkeleyGW, Copyright (c) 2011, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

Work done for this code was performed in the National University of Singapore (NUS). We ask that users of this code should cite the following reference:

> [K. Noori, N. L. Q. Cheng, F. Xuan, S. Y. Quek, 2D Materials **6**, 035036 (2019)](https://doi.org/10.1088/2053-1583/ab1e06)

If the BerkeleyGW code is used in conjunction with the wrrp.x code, please also cite the appropriate references for the BerkeleyGW code.

Authors: Keian Noori (c2dkn@nus.edu.sg), Nicholas Lin Quan Cheng (e0227689@u.nus.edu), Xuan Fengyuan (c2dxf@nus.edu.sg)

Contact: Quek Su Ying (phyqsy@nus.edu.sg)

## Features
- Real-space representation of the screening potential, inverse dielectric matrices, polarizability matrices in a user-specified supercell
- OpenMP parallelization

## Usage

### Compilation
wrrp.x can be compiled with both `ifort` and `gfortran`. The following dependencies are required:
- Intel MKL *or*
- fftw2 and LAPACK

To compile the complex version of the code using ifort and Intel MKL:
```
./cplx.sh
```
Note that you may need to modify this script to suit your machine's environment.

### Patching BerkeleyGW v1.2.0
In order to properly use the symmetry of the crystal (i.e. in order to use BerkeleyGW epsmat files generated in a reduced Brillouin zone), wrrp.x requires the user to extract symmetry information from the BerkeleyGW code. This ca be accomplished by using a modified version of BerkeleyGW v1.2.0, which writes two files, `isortg` and `gmapdata`. These files are written by the modified epsilon.x and sigma.x routines, respectively. We provide a patch file, `BerkeleyGW-1.2.0_qmap.patch`, in the `BerkeleyGW` directory, that applies the changes to BerkeleyGW v1.2.0 needed to write `isortg` and `gmapdata`. Upon entering the Berkeley v1.2.0 directory, the patch can be applied using the standard Unix 'patch' command:
```
cd /path/to/BGW1.2.0
patch -p0 < /path/to/patch/file/BerkeleyGW-1.2.0_qmap.patch
```
BerkeleyGW can then be compiled as usual.

### Running the code
#### Input files
- Input file `wrrp.inp`. See `wrrp.inp` in this directory for input parameters.
- q-point information `wrrp.qibz` and `wrrp.qfbz`. The input is extracted from the BerkeleyGW `kgrid.x` log file.
- Inverse dielectric matrices `epsmat` and `eps0mat` generated using BerkeleyGW.
-- Note that these **must** be in binary format (if compiled with HDF5 support, BerkeleyGW can force binary format output by using the `dont_use_hdf5` flag).
-- Note that we strongly recommend the use of BerkeleyGW v2.0 for the generation of `epsmat` files. The modified BerkeleyGW v1.2.0 should only be used to write `isortg`, `gmapdata`, and, if needed, `wrrp.wcoul0`.
- (If using Monte Carlo averaging to treat q -> 0 limit in W and V) `wrrp.wcoul0` (Refer to 1. in Notes)
- (If using nonuniform neck subsampling (NNS)) `subweights.dat` generated using BerkeleyGW
- (If using symmetry operations) `isortg` and `gmapdata`

#### Output files
- `data_rrp.bin` contains the final (6D) result in real space for further post-processing.
- (If selected in `wrrp.f90`) `wrrp.real|cplx|abs.out` provides the planar-averaged or 2D contour plot of the real, complex or absolute value of the output data.

#### Parallelization using OpenMP
- the code supports OpenMP parallelization via the `OMP_NUM_THREADS` variable
- on multicore machines, we recommend that users select a single MPI process and set `OMP_NUM_THREADS` according to their needs

### Notes
1. **Generation of `wrrp.wcoul0`** requires the user to explicitly write out the Monte Carlo averaged value of the q -> 0 limit in W and V during the computation of the self-energies at the Sigma step in BerkeleyGW. Refer to Computer Physics Communications 183 (2012) 1269–1289 for further discussion. The provided `BerkeleyGW-1.2.0_qmap.patch` enables Sigma to explicitly output `wrrp.wcoul0`. The patch works with BerkeleyGW/1.2.0.

2. **Symmetry operations.** The computation of the inverse dielectric matrices exploits symmetry operations to enable calculations within the irreducible Brillouin zone. In order to correctly include the symmetry operations in calculations using wrrp.x, the user is required to compute `isortg` and `gmapdata` at the Epsilon and Sigma step in BerkeleyGW, respectively. We provide a patch in `irrbz.patch` which enables the user to generate `isortg` and `gmapdata` without computing the inverse dielectric matrix and the self-energies respectively. The patch works with BerkeleyGW/1.2.0. Note that if the files `isortg` and `gmapdata` are not provided, the inverse dielectric matrix must be computed in the full Brillouin zone to give correct results.

### Known Issues
1. When the Monte-Carlo averaging scheme is used, the head element of W (written to `wrrp.wcoul0`) results in a planar-averaged Vscr(r,r) that does not tend toward 0 away from the slab. While the overall shape of the curve is correct, the offset is not. This issue is avoided if the NNS scheme is used. We recommend the use of NNS wherever possible.

2. In systems without symmetry (i.e. where the full BZ is used to calculate the wavefunctions), the `gmapdata` file might be incorrect, causing wrrp.x to crash.

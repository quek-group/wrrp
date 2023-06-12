# Autocorrelation
To understand the role of 2D materials in reducing charge
fluctuations, in `autocorr.f90` we compute the autocorrelation 
function C(r) of the screened Coulomb potential and bare 
Coulomb potential based on Eq. (16) in 2D Mater. 6 (2019) 035036.

Input data requires W(r,r') as computed from wrrp.x. 

Example for input file in autocorr.inp

# Wavefunction interpolation
To go beyond a point charge model to compute the effect of Vscr(r,r)
on a molecule, in `wfcint.f90` we use Bader charges centered at 
each atom of the molecule to improve our estimate of the HOMO-LUMO
gap of a molecule adsorbed on a 2D material.

Input data requires Vscr(r,r') as computed from wrrp.x. 

Example for input file in wfcint.inp

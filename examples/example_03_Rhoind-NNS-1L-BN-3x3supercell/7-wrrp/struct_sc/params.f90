!=========================================================================
! This module contains definitions for repeated parameters and constants
!=========================================================================

module params_m
  implicit none

  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: DPC = kind((1.0d0,1.0d0))
  real(DP), parameter :: PI_D = 3.1415926535897932384626433832795_dp
  real(DP), parameter :: E_D = 2.7182818284590452353602874713526_dp
  real(DP), parameter :: TOL_ZERO = 1.0d-12
  real(DP), parameter :: TOL_Small = 1.0d-6
  real(DP), parameter :: TOL_QG_ZERO = 0.003464102 ! corresponds to sqrt(0.002^2 + 0.002^2 + 0.002^2) 
  complex(DPC), parameter :: I_D = (0.0d0,1.0d0)
end module params_m

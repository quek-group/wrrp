!=========================================================================
! Module v
! Contains the routines necessary to calculate the Coulomb potential.
! 
!----------------------------
! version 0.1:  17.05.2016      
! KN: G space implementation
! version 0.2:  29.09.2016      
! KN: Corrected q + G addition
!=========================================================================

#include "f_defs.h"

module v_m

  use params_m
  implicit none

  public :: vcoul_g

contains

!------------------------------
! Calculate the Coulomb interaction v(G)=8*pi/|G|^2 (Ry)
! In the case of WGG', we need v(q+G').
  subroutine vcoul_g(gp,q,b,blat,vcoul)
    integer, intent(in) :: gp(3)
    real(DP), intent(in) :: q(3), b(3,3), blat ! G', q, b vectors, 2pi/alat
    SCALAR, intent(out) :: vcoul
    real(DP) :: qcart(3), gcart(3), vec(3), norm2

! Note that G vectors are in blat (2pi/alat) units and q vectors are in crystal
! coordinates.

    ! Convert q to Cartesian coordinates in blat units
    qcart = q(1)*b(:,1) + q(2)*b(:,2) + q(3)*b(:,3)

    ! Convert G' to Cartesian coordinates
    gcart = gp(1)*b(:,1) + gp(2)*b(:,2) + gp(3)*b(:,3)

    ! Compute sum of G' and q in blat units
    vec = 0.d0
    vec(1) = gcart(1)+qcart(1)
    vec(2) = gcart(2)+qcart(2)
    vec(3) = gcart(3)+qcart(3)

!!!!!LOG
     write(85,'(A,3f10.6)') "Bare Coulomb potential (vcoul) for q-point:",q
     write(85,'(A,3i5.1)') "and G' vector:",gp
     write(85,'(A,3f10.6)') "q in Cartesian 2pi/a units is:",qcart
     write(85,'(A,3f10.6)') "G' in Cartesian 2pi/a units is:",gcart
     write(85,'(A,3f10.6)') "The q+G' vector in Cartesian 2pi/a units is:", vec
!!!!!LOG

    ! Convert to 1/au units
    vec(1) = vec(1)*blat
    vec(2) = vec(2)*blat
    vec(3) = vec(3)*blat
    norm2 = vec(1)**2 + vec(2)**2 + vec(3)**2 

!!!!!LOG
     write(85,'(A,3f10.6)') "The q+G' vector in Cartesian 1/au units is:", vec
!!!!!LOG

! Calculate Coulomb interaction. 
! Note that we use 8*pi insteady of 4*pi since we want Ry and not Ha

!KN: note that for now we don't use Coulomb trucation
!KN: note that for now we set v[(q+G')=0] = 0 
    if (norm2 .eq. 0) then
      vcoul = 0
    else 
      vcoul = (8*PI_D)/norm2
    endif

!!!!!LOG
    write(85,'(A,f20.6)') "v(q+G') in Ry is", vcoul
    write(85,*) '--------------------'
!!!!!LOG

    return
  end subroutine vcoul_g

end module v_m

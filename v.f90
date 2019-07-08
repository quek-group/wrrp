!=========================================================================
! Module v
! Contains the routines necessary to calculate the bare Coulomb potential
!
! This code was written at the National University of Singapore (NUS)
! v0.1 K Noori
!=========================================================================

#include "f_defs.h"

module v_m

  use params_m
  implicit none

  public :: vcoul_g
  public :: vcoul_ggp

contains

  !------------------------------
  ! Calculate the asymmetric  Coulomb interaction
  ! v(G)=8*pi/|G|^2 (Ry)
  ! In the case of WGG', we need v(q+G').
  subroutine vcoul_g(use_slab,gp,q,b,blat,vcoul0,vcoul)
    logical, intent(in) :: use_slab                     ! 2D slab truncation
    integer, intent(in) :: gp(3)
    real(DP), intent(in) :: q(3), b(3,3), blat, vcoul0  ! G', q, b vectors, 2pi/alat, V(q->0,G=0) from BGW
    SCALAR, intent(out) :: vcoul
    real(DP) :: qcart(3), gcart(3), qg(3), norm_cart, norm_sq, qgxy(3), qgz(3), kxy, kz, zc, tol_head

    ! Note that G vectors are in blat (2pi/alat) units and q vectors are in crystal
    ! coordinates.

    ! Convert q to Cartesian coordinates in blat units
    qcart = q(1)*b(:,1) + q(2)*b(:,2) + q(3)*b(:,3)
    ! Convert G' to Cartesian coordinates
    gcart = gp(1)*b(:,1) + gp(2)*b(:,2) + gp(3)*b(:,3)
    ! Compute sum of G' and q in blat units
    qg = 0.d0
    qg = gcart + qcart
    ! Magnitude of q+G in blat units
    norm_cart = norm2(qg)

    ! Set the threshold at which we consider |q+G| = 0, i.e. head element
    ! This sets the threshold below at or below which we use vcoul0
    !! Default is TOL_QG_ZERO, which corresponds to a maximum q0 of
    !! (0.002,0.002,0.002) in blat units
    !!tol_head = TOL_QG_ZERO
    ! KN: We're now directly passing a q = 0 vector when in the q0 loop (aside
    !     from truncated metals).
    !     As a result we want use vcoul0 when (q+G)=0 (i.e. when |q+G|
    !     < TOL_ZERO) in order to avoid a divergence.
    tol_head = TOL_ZERO
    !    tol_head = TOL_QG_ZERO

!!!!!LOG
    ! Nic: I find this logging to be time-consuming...
    !      maybe some verbosity flag to enable this?
    !    write(85,'(a50,3f10.6)') "Bare Coulomb potential (vcoul) for q-point:",q
    !    write(85,'(a50,3i10.1)') "and G' vector:",gp
    !    write(85,'(a50,3f10.6)') "q in Cartesian blat units is:",qcart
    !    write(85,'(a50,3f10.6)') "G' in Cartesian blat units is:",gcart
    !    write(85,'(a50,3f10.6)') "The q+G' vector in Cartesian blat units is:", qg
!!!!!LOG

    ! Convert to 1/au units
    qg = qg * blat
    ! Compute squared norm of (q+G')
    norm_sq = norm2(qg)**2

!!!!!LOG
    ! Nic: time-consuming logging?
    !    write(85,'(a50,3f10.6)') "The q+G' vector in Cartesian 1/bohr units is:", qg
!!!!!LOG

    ! Calculate Coulomb interaction
    ! Note that we use 8*pi insteady of 4*pi since we want Ry and not Ha

    if (norm_cart .le. tol_head) then
       ! This is the head element (|q+G| = 0) so use vcoul0
       write(6,*) "Replacing head element of V with vcoul0..."
       vcoul = vcoul0
    elseif (use_slab) then
       !if (use_slab) then
       ! Slab Coulomb truncation is enabled
       ! Use Ismail-Beigi's analytical expression, as implement in Berkeley GW
       ! non-periodic direction MUST be along z (i.e. b(:,3))
       qgxy(1:2) = qg(1:2)
       qgxy(3) = 0.d0
       kxy = sqrt(dot_product(qgxy,qgxy))
       kz = qg(3)
       zc = (PI_D)/(blat*b(3,3))
       vcoul = (8.0d0*PI_D)/norm_sq
       vcoul = vcoul*(1.0d0-exp(-kxy*zc)*cos(kz*zc))
    else
       ! No Coulomb truncation
       vcoul = (8*PI_D)/norm_sq
    endif

!!!!!LOG
    ! Nic: time-consuming logging
#ifdef CPLX
    !    write(85,'(a50,2f14.8)') "v(q+G') in Ry is:", vcoul
#else
    !    write(85,'(a50,f14.8)') "v(q+G') in Ry is:", vcoul
#endif
    !    write(85,*) repeat("=",79)
!!!!!LOG

    return
  end subroutine vcoul_g

  !------------------------------
  ! Calculate the symmetric  Coulomb interaction
  ! v(G,G')=8*pi/(|G||G'|) (Ry)
  ! In the case of WGG', we need v(q+G,q+G').
  subroutine vcoul_ggp(g,gp,q,b,blat,vcoul)
    integer, intent(in) :: g(3), gp(3)
    real(DP), intent(in) :: q(3), b(3,3), blat ! G', q, b vectors, 2pi/alat
    SCALAR, intent(out) :: vcoul
    real(DP) :: qcart(3), gcart(3), gpcart(3), qg1(3), qg2(3)

    ! Note that G vectors are in blat (2pi/alat) units and q vectors are in crystal
    ! coordinates.

    ! Convert q to Cartesian coordinates in blat units
    qcart = q(1)*b(:,1) + q(2)*b(:,2) + q(3)*b(:,3)

    ! Convert G to Cartesian coordinates
    gcart = g(1)*b(:,1) + g(2)*b(:,2) + g(3)*b(:,3)

    ! Convert G' to Cartesian coordinates
    gpcart = gp(1)*b(:,1) + gp(2)*b(:,2) + gp(3)*b(:,3)

    ! Compute sum of G+q and G'+q in blat units
    qg1 = 0.d0
    qg1 = qcart + gcart
    qg2 = 0.d0
    qg2 = qcart + gpcart

    ! Convert to 1/au units
    qg1 = qg1 * blat
    qg2 = qg2 * blat

    ! Calculate Coulomb interaction.
    ! Note that we use 8*pi insteady of 4*pi since we want Ry and not Ha

    !KN: note that for now we don't use Coulomb trucation
    !KN: note that for now we set v = 0 if q+G or q+G' are 0
    if (norm2(qg1) .eq. 0 .or. norm2(qg2) .eq. 0) then
       vcoul = 0
    else
       vcoul = (8*PI_D)/(norm2(qg1)*norm2(qg2))
    endif

    return
  end subroutine vcoul_ggp

end module v_m

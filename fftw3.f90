!=========================================================================
! Module fftw
! Contains the routines necessary to setup and perform FFTW calculations. 
! Adapted from fftw.f90 of Berkeley GW.
! 
!----------------------------
! version 0.1:  18.04.2016      
! KN: Initial test
!----------------------------
! version 0.2:  03.05.2016
! KN: Included definition of SCALAR via f_defs.h
!----------------------------
! version 0.3:  05.07.2016
! KN: Calculation of real space grid.
! Nic: Using fftw3 if want to exploit threading.
!=========================================================================

#include "f_defs.h"

module fftw_m

  use params_m
  use supercell_m
  use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f'
  integer*8, private :: fft_plan = 0
  integer, private :: ifirst = 0
  integer, private :: num_fft_threads = 0
  integer :: iret
  private :: check_fft_size,   &
             gvec_to_fft_index 

  public  :: setup_fft_sizes,  &
             put_into_fftbox,  &
             get_from_fftbox,  &
             get_real_grid,    &
             do_fft
!----------------------------
! Interfaces to FFTW routines.
! Taken from Berkeley GW.

  integer, private :: nfftold(6) = 0 ! Old FFT box dimensions


contains

!------------------------------
! Calculate the ideal FFT box size based on max dimension of FFT grid.
! scaling value set in order to normalize data since FFTW doesn't scale
  subroutine setup_fft_sizes(nmax,nfft,scale)
    integer, intent(in) :: nmax(3)
    integer, intent(out) :: nfft(6)
    real(DP), intent(out) :: scale
    
    integer, parameter :: nfac = 3
    integer :: i

    ! Set ideal box size of FFT grid for G
    do i = 1,3
      nfft(i) = nmax(i)
      do while (.not. check_fft_size(nfft(i),nfac))
        nfft(i) = nfft(i) + 1
      enddo
    enddo
    ! Ideal FFT dimensions for G' are the same as for G     
    nfft(4) = nfft(1)
    nfft(5) = nfft(2)
    nfft(6) = nfft(3)

! KN: calculating scale here may be unnecessary due to maxnfft in wrrp.f90    
    scale = 1.0d0/product(nfft(1:6))
    return
  end subroutine setup_fft_sizes

!------------------------------  
  logical function check_fft_size(nfft,nfac)
    integer, intent(in) :: nfft, nfac
 
    integer, parameter :: maxfac = 6
    integer, parameter :: fac(maxfac) = (/ 2, 3, 5, 7, 11, 13 /) 
    integer :: remainder, product, ifac, ipow, maxpow, pow(maxfac)

! KN: Note that nfft and nfac must be .ge. 1 and nfac must be .le. maxfac.
!      No explicit check is made here...at least not yet.

    remainder = nfft 
    do ifac = 1, maxfac
      pow(ifac) = 0 
    enddo
    
    do ifac = 1, nfac
      maxpow = int(log(dble(remainder)) / log(dble(fac(ifac)))) + 1
      do ipow = 1, maxpow
        if (mod(remainder, fac(ifac)) .eq. 0) then 
          remainder = remainder / fac(ifac)
          pow(ifac) = pow(ifac) + 1
        endif
      enddo
    enddo
    
    product = remainder
    do ifac = 1, nfac
      do ipow = 1, pow(ifac)
        product = product * fac(ifac)
      enddo 
    enddo 

! KN: Should add here error handling
!      if (product .ne. nfft) then DIE: factorization failed

    check_FFT_size = remainder .eq. 1 .and. pow(5) .le. 1 .and. pow(6) .le. 1
  end function check_fft_size

!------------------------------
! This subroutine takes the KE-sorted G vectors G(1:3) and G'(1:3) and the 6D
! FFT box size nfft(1:6) and finds the points in the box corresponding to those
! vectors.
! This routine modifies the corresponding Berkeley GW routine in order to handle
! a 6D FFT.
  subroutine gvec_to_fft_index(g,gp,ig,igp,nfft)
    integer, intent(in) :: g(3), gp(3), nfft(6)
    integer, intent(out) :: ig(3)
    integer, intent(out) :: igp(3)

    ig(1:3) = g(1:3) + 1
    igp(1:3) = gp(1:3) + 1

    if (g(1) < 0) ig(1) = nfft(1) + ig(1)
    if (g(2) < 0) ig(2) = nfft(2) + ig(2)
    if (g(3) < 0) ig(3) = nfft(3) + ig(3)
    if (gp(1) < 0) igp(1) = nfft(4) + igp(1)
    if (gp(2) < 0) igp(2) = nfft(5) + igp(2)
    if (gp(3) < 0) igp(3) = nfft(6) + igp(3)

    return
  end subroutine gvec_to_fft_index

!------------------------------
! This routine takes data(1:ndata) and puts it into the 6D FFT box fftbox(:,:,:,:,:,:).
! This routine modifies the corresponding Berkekely GW routine in order to
! handle a 6D FFT.
! 
!   ndata -- number of data items in data(:)
!   data -- the data set, real or complex, depending on ifdef CPLX
!   ngused -- number of G vectors used in epsilon calculation (i.e. below eps
!             cutoff)
!   ngtot -- total number of G vectors in glist
!   glist -- a master list of G vectors
!   isortg(1:ng) -- correspondence between KE-sorted G vectors - used to store epsilon - and                    the unsorted G vectors in glist
!   fftbox(:,:,:,:,:,:) -- 6D complex FFT box where the data is put
!   nfft(1:6) -- sizes of FFT box: Nx,Ny,Nz,Npx,Npy,Npz
  subroutine put_into_fftbox(data, ngused, ngtot, glist, fftbox, nfft, isortg)
    SCALAR, intent(in) :: data(:) ! (ngused^2)
    integer, intent(in) :: ngused 
    integer, intent(in) :: ngtot
    integer, intent(in) :: glist(:,:) ! (3, ngtot)
    integer, optional, intent(in) :: isortg(:) ! (ngtot)
    integer, intent(in) :: nfft(:) ! (6)
    complex(DPC), intent(out) :: fftbox(:,:,:,:,:,:) ! (nfft(1), nfft(2), ... ,nfft(6))

    integer :: i, j, k, idata, iboxg(3), iboxgp(3), g(3), gp2(3)
    integer :: ig, igp

    ! Zero FFT box and put data into it
    fftbox(:,:,:,:,:,:) = (0.0d0,0.0d0)
    idata = 1 ! from 1 to ngused^2

!$omp parallel do private(g,gp2,iboxg,iboxgp)
    do i=1,ngused
      do j=1,ngused
        ! In the case of calctype = 4 (chi) we don't use a KE sorting of the G
        ! vectors, therefore isortg will not be present and we don't map the G
        ! vector indices
        if (present(isortg)) then
          g = glist(:,isortg(i))
          gp2 = -glist(:,isortg(j))
        else
          g = glist(:,i)
          gp2 = -glist(:,j)
        endif
        call gvec_to_fft_index(g,gp2,iboxg,iboxgp,nfft)
        fftbox(iboxg(1),iboxg(2),iboxg(3),iboxgp(1),iboxgp(2),iboxgp(3)) = data(idata)
        idata = idata + 1
      enddo
    enddo
!$omp end parallel do
    return
  end subroutine put_into_fftbox

!------------------------------
! Does the inverse of put_into_fftbox routine. Takes data in fftbox(:,:,:,:,:,:)
! and puts it back into the 1D data(1:ngused^2) array. The routine is slightly
! modified from the corresponding Berkeley GW routine to accomodate the 6D FFT
! box. Note that since FFTW does not scale data after an FFT, the data in the
! box is multiplied by 'scale' before being stored in data().
  subroutine get_from_fftbox(data, ngused, ngtot, glist, gindex, fftbox, nfft, scale)
    SCALAR, intent(out) :: data(:) ! (ngused^2)
    integer, intent(in) :: ngused
    integer, intent(in) :: ngtot
    integer, intent(in) :: glist(:,:) ! (3, ngtot)
    integer, intent(in) :: gindex(:) ! (ngtot)
    integer, intent(in) :: nfft(:) ! (6)
    complex(DPC), intent(in) :: fftbox(:,:,:,:,:,:) ! (nfft(1), nfft(2), ... ,nfft(6))
    real(DP), intent(in) :: scale

    integer :: i, j, idata, iboxg(3), iboxgp(3)

    ! Zero data vector and put back data from FFT box
    data(:) = 0.0
    idata = 1 ! from 1 to ngused^2
    do i=1,ngused
      do j=1,ngused
        call gvec_to_fft_index(glist(:,gindex(i)),glist(:,gindex(j)),iboxg,iboxgp,nfft)     
        data(idata) = fftbox(iboxg(1),iboxg(2),iboxg(3),iboxgp(1),iboxgp(2),iboxgp(3))*scale
        idata = idata + 1
      enddo
    enddo

    return
  end subroutine get_from_fftbox

!------------------------------
! Determines the real space grid (in Cartesian coordinates) from the real space lattice 
! vectors and the size of the FFT box, nfft.
subroutine get_real_grid(alat,a,nfft,rcart,rvec2real,sc_size,sc2uc_idx, &
                           scvec2real,sc_rcart,sc_a)
    integer, intent(in) :: nfft(6)
    real(DP), intent(in) :: alat, a(3,3) ! in units of alat
    real(DP), intent(out) :: rcart(:,:) ! (nfft(1)*nfft(2)*nfft(3),3)
    integer, intent(out) :: rvec2real(:,:) ! maps r vectors to real box index
    integer :: i, j, k, ir, allocerr, nr
!Nic
    integer, intent(in) :: sc_size(3)
    integer, intent(out) :: sc2uc_idx(:), scvec2real(:,:)
    real(DP), intent(out) :: sc_rcart(:,:), sc_a(3,3)
    integer :: nmax
!
    ! real space point intervals (deltas) for each lattice vector
    real(DP) :: dela1(3), dela2(3), dela3(3)
!    real(DP), allocatable :: r(:,:) ! (3,nfft(1)*nfft(2)*nfft(3))
!    integer, allocatable :: r2real ! (3,nr) - temp. holder for rvec2real

    nr = nfft(1)*nfft(2)*nfft(3) ! number of real space points
    nmax = nr*sc_size(1)*sc_size(2)*sc_size(3) ! real space pts in supercell
!    allocate(r (3,nr), STAT=allocerr)
!    allocate(r2real (3,nr), STAT=allocerr)
    dela1 = a(:,1)/nfft(1)
    dela2 = a(:,2)/nfft(2)
    dela3 = a(:,3)/nfft(3)

    ! determine all real space points in a() basis in units of alat
    ir = 1
    do i = 1,nfft(1)
      do j = 1,nfft(2)
        do k = 1,nfft(3)
          ! The first FFTW output is 0
          rcart(:,ir) = (i-1)*dela1 + (j-1)*dela2 + (k-1)*dela3
          rvec2real(1,ir) = i
          rvec2real(2,ir) = j
          rvec2real(3,ir) = k
          ir = ir + 1
        enddo
      enddo
    enddo
    
    call get_sc_grid(a,nfft,rcart,sc_size,sc_a,sc2uc_idx,scvec2real,sc_rcart)
!!!!!LOG
    write(85,*) "----------------------------------------------"
    write(85,*) "The real space grid deltas are:"
    write(85,*) (dela1(i),i=1,3)
    write(85,*) (dela2(i),i=1,3)
    write(85,*) (dela3(i),i=1,3)
    write(85,*) "----------------------------------------------"
    write(85,*) "The real space Cartesian points in alat are:"
    write(85,*) ((rcart(i,j),i=1,3),j=1,nr)
    write(85,*) "----------------------------------------------"
    write(85,*) "The real space Cartesian points in alat are:"
    do j=1,nmax
    write(85,*) (sc_rcart(i,j),i=1,3), (scvec2real(i,j),i=1,3)
    enddo
!!!!!LOG
    
!    ! transform basis from a() to Cartesian
!    do ir = 1,nr
!      rcart(1,ir) = r(1,ir)*a(1,1) + r(2,ir)*a(1,2) + r(3,ir)*a(1,3)
!      rcart(2,ir) = r(1,ir)*a(2,1) + r(2,ir)*a(2,2) + r(3,ir)*a(2,3)
!      rcart(3,ir) = r(1,ir)*a(3,1) + r(2,ir)*a(3,2) + r(3,ir)*a(3,3)
!    enddo
!
!!!!!!LOG
!    write(85,*) "----------------------------------------------"
!    write(85,*) "The real space Cartesian points in alat are:"
!    write(85,*) ((rcart(i,j),i=1,3),j=1,nr)
!!!!!!LOG

    ! transform units from alat to au
!    rcart = r
    rcart = alat * rcart
    sc_rcart = alat * sc_rcart
!!!!!LOG
!    write(85,*) "----------------------------------------------"
!    write(85,*) "The real space Cartesian points in bohr are:"
!    write(85,*) ((rcart(i,j),i=1,3),j=1,nr)
    write(85,*) "----------------------------------------------"
    write(85,*) "The real space Cartesian points in bohr are:"
    write(85,*) ((sc_rcart(i,j),i=1,3),j=1,nmax)
!!!!!LOG
    
    return
  end subroutine get_real_grid

!------------------------------
  subroutine do_fft(fftbox, nfft, sign)
    complex(DPC), intent (in) :: fftbox(:,:,:,:,:,:)
    integer, intent (in) :: nfft(6), sign

! KN: Variable sign must be either -1 (forward FFT) or 1 (backward FFT).
!     No explicit check is made at the moment.

! Create the FFTW plans

! KN: Here there should be a check to see if there already exists a plan for a
!     given nfft(). If so it shouldn't be recreated. Check via nfftold variable.

    if (ifirst .eq. 0) then
      call dfftw_init_threads(iret)
      ifirst = 1
!$OMP PARALLEL
      num_fft_threads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
      call dfftw_plan_with_nthreads(num_fft_threads)
    endif

    if (sign == 1) then
      call dfftw_plan_dft(fft_plan,6,nfft,fftbox,fftbox,FFTW_BACKWARD,FFTW_ESTIMATE)
    else if (sign == -1) then
      call dfftw_plan_dft(fft_plan,6,nfft,fftbox,fftbox,FFTW_FORWARD,FFTW_ESTIMATE)
    endif
    call dfftw_execute_dft(fft_plan,fftbox,fftbox)
    call dfftw_destroy_plan(fft_plan)
    nfftold(:) = -1
    return
  end subroutine do_fft

end module fftw_m

!=========================================================================
! Module input
! Contains routines relating to the input file for the program.
! Adapted from Berkeley GW inread subroutine.

! This code was written at the National University of Singapore (NUS)
! v0.1 K Noori
! v0.2 K Noori, N Cheng
!=========================================================================

#include "f_defs.h"

module input_m

  use params_m
  use fftw_m
  implicit none

  public  :: read_input
  public  :: read_q_points
  public  :: read_qs_weights
  public  :: read_wcoul0
  public  :: binaryfile_init
  public  :: binaryfile_setup
  public  :: binaryfile_readmatrix

contains

  !------------------------------
  ! Read the input file 'wrrp.inp'
  subroutine read_input(a,b,alat,blat,vol,use_q0,use_ss,use_wavg,use_slab,calctype,calcname,sc_size)
    character, intent(out) :: calcname*12
    integer, intent(out) :: calctype ! calculation type fo wrrp.x
    real(DP), intent(out) :: a(3,3), b(3,3), alat, blat, vol
    logical, intent(out) :: use_q0, use_ss, use_wavg, use_slab ! q0 point indicator, W head averaging, slab truncation
    integer :: ioerr, i, j
    !Nic
    integer, intent(out) :: sc_size(3) ! supercell size

    ! KN: Read the input file, hardcoded as 'wrrp.inp'
    !     File is assumed to exist; no error checking done yet.
    ! Nic:status='old' means that if file doesn't exist
    !     the program will just terminate violently
    open(unit=77,file='wrrp.inp',form='formatted',status='old')

    ! Set default values KN: Currently values in wrrp.inp must be set
    ! in correct order This should be changed to a keyword system
    calctype = 0        ! calculation type
    use_slab = .false.      ! use 2D (slab) Coulomb truncation for V
    ! KN: FIX THIS: use_wavg should be required (i.e. 'wrrp.wcoul0' must be present)
    ! unless use_ss = .true., in which case it should not be used
    use_wavg = .false.      ! use MC averaging for W and V as q-->0
    use_q0 = .true.     ! indicates the presence of a BGW eps0mat file
    ! KN: FIX THIS: the use of subsampling only makes sense when computing W or Vscr.
    ! setting use_ss = .true. for calctype other than 0 or 1 should cause
    ! the program to die gracefully
    use_ss = .false.    ! use BGW subsampling scheme for q-->0


    ! Read calculation type
    ! 0: screened coulomb: W (default)
    ! 1: screening potential: Vscr = W - V
    ! 2: bare coulomb: V
    ! 3: inverse dialectric function: epsilon^-1
    ! 4: polarizability: chi
    ! 5: induced charge: rhoind
    read(77,*,iostat=ioerr) calctype
    if (calctype .eq. 0) then
       calcname = "W(r,r')"
    elseif (calctype .eq. 1) then
       calcname = "Vscr(r,r')"
    elseif (calctype .eq. 2) then
       calcname = "V(r)"
    elseif (calctype .eq. 3) then
       calcname = "epsinv(r,r')"
    elseif (calctype .eq. 5) then
       calcname = "rhoind(r,r')"
    else
       calcname = "chi(r,r')"
    endif

    ! Read slab truncation boolean
    ! .false. (default): no 2D slab truncation
    ! .true.: enable 2D slab truncation
    read(77,*,iostat=ioerr) use_slab

    ! Read use_wavg boolean
    ! .true. : replace V(q->0,G=0) and W(q->0,G=0) with vcoul0 and wcoul0,
    ! respectively, from Berkeley GW via wrrp.wcoul0 input file
    read(77,*,iostat=ioerr) use_wavg

    ! Read use_q0 boolean
    ! .true. (default): read q-> 0 eps0mat/chi0mat file
    ! .false.: do not read eps0mat/chi0mat file
    read(77,*,iostat=ioerr) use_q0

    ! Read use_ss boolean
    ! .false. (default): do not read 'subweights.dat' file
    ! .true.: read BGW-generated 'subweights.dat' file
    read(77,*,iostat=ioerr) use_ss

    ! Read lattice vectors a
    do i=1,3
       read(77,*,iostat=ioerr) (a(j,i),j=1,3)
    enddo

    ! Read supercell dimensions
    read(77,*,iostat=ioerr) (sc_size(i),i=1,3)

    ! Direct lattice vector magnitude alat
    alat = sqrt(a(1,1)**2 + a(2,1)**2 + a(3,1)**2) ! units of au
    ! Reciprocal lattice vector magnitude
    blat = (2.0d0 * PI_D) / alat
    ! Normalize lattice vectors to units of alat
    a(:,:) = a(:,:) / alat
    ! Cell volume in alat^3
    vol =  a(1,1) * (a(2,2) * a(3,3) - a(2,3) * a(3,2)) - &
         a(2,1) * (a(1,2) * a(3,3) - a(1,3) * a(3,2)) + &
         a(3,1) * (a(1,2) * a(2,3) - a(1,3) * a(2,2))

    ! Determine reciprocal lattice vectors b in units of 2pi/alat
    b(1,1) = (a(2,2) * a(3,3) - a(3,2) * a(2,3)) / vol
    b(2,1) = (a(3,2) * a(1,3) - a(1,2) * a(3,3)) / vol
    b(3,1) = (a(1,2) * a(2,3) - a(2,2) * a(1,3)) / vol
    b(1,2) = (a(2,3) * a(3,1) - a(3,3) * a(2,1)) / vol
    b(2,2) = (a(3,3) * a(1,1) - a(1,3) * a(3,1)) / vol
    b(3,2) = (a(1,3) * a(2,1) - a(2,3) * a(1,1)) / vol
    b(1,3) = (a(2,1) * a(3,2) - a(3,1) * a(2,2)) / vol
    b(2,3) = (a(3,1) * a(1,2) - a(1,1) * a(3,2)) / vol
    b(3,3) = (a(1,1) * a(2,2) - a(2,1) * a(1,2)) / vol

    close(77)

  end subroutine read_input

  !------------------------------
  ! Read the wrrp.qibz and wrrp.qfbz input files in order to get the list of q-points in
  ! the full and irreducible BZ, along with the mapping between the two
  subroutine read_q_points (qfbz,nqfbz,qmapf,qmapi)
    integer :: i, nqibz, ioerr, allocerr
    real(DP) :: tmp
    character :: dummy
    real(DP), allocatable :: wq(:) ! weights of the q-pts in the irreducible BZ; not currently needed
    integer, intent(out) :: nqfbz ! number of q in full BZ
    integer, allocatable, intent(out) :: qmapi(:) ! mapping of index of q-pt in iBZ to its index in fBZ
    integer, allocatable, intent(out) :: qmapf(:) ! mapping of a folded fBZ q-pt to corresponding unfolded one
    real(DP), allocatable, intent(out) :: qfbz(:,:) ! q-pts in full BZ, weights in the irreducible BZ

    open(unit=36,file='wrrp.qfbz',form='formatted',status='old',action='read',iostat=ioerr)
    open(unit=37,file='wrrp.qibz',form='formatted',status='old',action='read',iostat=ioerr)

    read(36,*) nqfbz ! Read total number of q points in full BZ
    allocate(qfbz (3,nqfbz), stat=allocerr)
    allocate(qmapf (nqfbz), stat=allocerr)
    do i = 1,nqfbz
       read(36,*,iostat=ioerr) tmp, qfbz(:,i), tmp, qmapf(i), dummy
    enddo

    read(37,*) nqibz ! Number of q points in irreducible BZ including q0, if present
    allocate(wq (nqibz), stat=allocerr)
    allocate(qmapi (nqibz), stat=allocerr)
    do i = 1,nqibz
       read(37,*,iostat=ioerr) tmp, tmp, tmp, tmp, wq(i), qmapi(i)
    enddo

    close(36)
    close(37)
  end subroutine read_q_points

  !------------------------------
  ! Read the input file 'subweights.dat' produced by the BGW subsample.x program
  ! in order to get the weights of the subsampled q-pts (qs) in the q=0 Voronoi cell.
  ! In 'subsample.dat' the first line lists the total number of qs. Each
  ! subsequent line gives the weight and magnitude of the corresponding qs
  subroutine read_qs_weights(nqs,wqs)
    integer, intent(out) :: nqs     ! number of subsampled (qs) q points
    integer :: iqs, allocerr, ioerr
    real(DP) :: val
    real(DP), allocatable, intent(out) :: wqs(:)     ! weights of the qs points

    open(unit=81,file='subweights.dat',form='formatted',status='old',iostat=ioerr)
    read(81,*,iostat=ioerr) nqs ! Read total number of qs points
    allocate(wqs (nqs),stat=allocerr)
    do iqs = 1,nqs
       read(81,*,iostat=ioerr) wqs(iqs)
       !if (ioerr .ne. 0) exit
       !iq = iq + 1
       !wq(iq) = val
    enddo
    close(81)
  end subroutine read_qs_weights

  !------------------------------
  ! Read the input file 'wrrp.wcoul0' produced by our modified BGW program.
  ! This file is needed when use_wavg is true and contains the screening type and
  ! head elements of bare and screened Coulomb (V & W) matrices determined in BGW
  ! by Monte Carlo averaging.
  subroutine read_wcoul0(iscreen,vcoul0,wcoul0)
    character :: dummy
    integer :: ioerr
    real(DP) :: wcoul0r, wcoul0i
    integer, intent(out) :: iscreen     ! BGW screening type
    real(DP), intent(out) :: vcoul0
    SCALAR, intent(out) :: wcoul0

    ! From BGW: iscreen:
    ! 0 = semiconductor
    ! 1 = graphene
    ! 2 = metal

    open(65,file='wrrp.wcoul0',form='formatted',status='old',iostat=ioerr)
    !      open(unit=121,file="slab.log",form='formatted',status='replace')
    read(65,*) dummy, iscreen ! BGW screening type
    read(65,*) !q0 vector
    read(65,*) !G0
    read(65,*) dummy, vcoul0
#ifdef CPLX
    read(65,*) dummy, wcoul0r, wcoul0i
    wcoul0 = cmplx(wcoul0r,wcoul0i)
#else
    read(65,*) dummy, wcoul0
#endif
    close(65)

  end subroutine read_wcoul0

  !------------------------------
  ! DEPRECATED -> REMOVE ME!
  ! Initial read of BerkeleyGW binary epsmat or chimat file in order
  ! to get number of q points and G vectors
  subroutine binaryfile_init(nq,ng,calctype)
    integer, intent(in) :: calctype
    integer, intent(out) :: ng, nq
    integer :: ioerr, nq0

    ! KN: Read epsmat/chimat file, code adapted directly from BerkeleyGW.
    !     This assumes freq_dep .eq. 0 (static/PPM).
    !     Should later account for freq_dep .eq. 2 (full freq) or 3 (GN PPM)
    if (calctype .ne. 4) then
       ! Read epsmat file
       open(unit=11,file='epsmat',form='unformatted',status='old',iostat=ioerr)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11) nq
       read(11) ng
       close(11)
    else
       ! Read chimat file
       open(unit=09,file='chimat',form='unformatted',status='old',iostat=ioerr)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09) ng
       read(09) nq, nq0
       close(09)
    endif

    ! In chimat nq includes q0 point, whereas in epsmat nq does not
    if (calctype .eq. 4) then
       nq = nq-nq0
    endif

  end subroutine binaryfile_init
  ! DEPRECATED -> REMOVE ME!

  !------------------------------
  ! Read BerkeleyGW binary epsmat or chimat file and extract
  ! q points, G vectors, and setup FFT mesh
  subroutine binaryfile_setup (nq,nqs,ng,kx,ky,kz,q,qs,q0,use_q0,use_ss,nfft,maxnfft,calctype)
    logical, intent(in) :: use_q0, use_ss
    integer, intent(in) :: calctype
    integer, intent(out) :: nfft(:), maxnfft(:), nq, nqs, ng
    integer, allocatable, intent(out) :: kx(:), ky(:), kz(:)
    real(DP), intent(out) :: q0(:)
    real(DP), allocatable, intent(out) :: q(:,:),qs(:,:)
    integer :: ioerr, allocerr, i, iq, ig, nq0, nmtx, maxdim, kmax(3), nmax(3)
    integer :: gridpts, grid(3), ngq ! dummy arguments
    real(DP) :: scale
    integer, allocatable :: isortg(:), isortgi(:) ! mapping between unordered and KE-ordered G vectors, and vice versa
    real(DP), allocatable :: q_dummy(:,:) ! dummy array used in allocation of q(:,:) for chimat files

    ! Read q points and G vectors from appropriate binary file
    if (calctype .ne. 4) then
       ! We are not looking at chi: read epsmat file

       if (use_q0) then
          ! Read q0 from eps0mat if use_q0 = .true.
          open(unit=12,file='eps0mat',form='unformatted',status='old',iostat=ioerr)
          if (use_ss) then
             ! We are using NNS subsampling scheme:
             ! read multiple qs-points from eps0mat
             !         open(unit=12,file='eps0mat',form='unformatted',status='old',iostat=ioerr)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12) nqs
             backspace(12)
             allocate (qs (3,nqs),stat=allocerr)
             read(12) nqs,((qs(i,iq),i=1,3),iq=1,nqs)
             ! Set q0 = qs(1)
             q0 = qs(:,1)
             !          read(12) ng
             !          backspace(12)
             !          allocate (kx (ng), stat=allocerr)
             !          allocate (ky (ng), stat=allocerr)
             !          allocate (kz (ng), stat=allocerr)
             !          read(12) ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)
             !         close(12)
          else
             ! We are NOT using the NNS subsampling scheme
             ! read a single q0 point from eps0mat
             !         open(unit=12,file='eps0mat',form='unformatted',status='old',iostat=ioerr)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12)
             read(12) nq,q0
             !         close(12)
          endif
          close(12)
       else
          ! Otherwise set q0 = (0,0,0)
          q0 = 0.d0
       endif

       ! Read main epsmat file in order to setup G vectors, remaning q-points,
       ! and FFT grid
       open(unit=11,file='epsmat',form='unformatted',status='old',iostat=ioerr)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11)
       read(11) nq
       backspace(11)
       allocate (q (3,nq),stat=allocerr)
       read(11) nq,((q(i,iq),i=1,3),iq=1,nq)
       read(11) ng
       backspace(11)
       allocate (kx (ng), stat=allocerr)
       allocate (ky (ng), stat=allocerr)
       allocate (kz (ng), stat=allocerr)
       read(11) ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)

    else
       ! We are looking at chi: read chimat file

       if (use_q0) then
          ! Read q0 from chi0mat if use_q0 = .true.
          ! Otherwise set q0 = (0,0,0)
          open(unit=10,file='chi0mat',form='unformatted',status='old',iostat=ioerr)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10)
          read(10) (q0(i),i=1,3)
          close(10)
       else
          q0 = 0.d0
       endif

       ! Read main chimat file in order to setup G vectors, remaning q-points,
       ! and FFT grid
       open(unit=09,file='chimat',form='unformatted',status='old',iostat=ioerr)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09)
       read(09) ng
       backspace(09)
       allocate (kx (ng), stat=allocerr)
       allocate (ky (ng), stat=allocerr)
       allocate (kz (ng), stat=allocerr)
       read(09) ng, gridpts, grid, (kx(ig),ky(ig),kz(ig),ig=1,ng)
       read(09) nq
       backspace(09)
       allocate (q_dummy (3,nq), stat=allocerr)
       read(09) nq, nq0, ((q_dummy(i,iq),i=1,3),iq=1,nq)

       ! nq in epsmat does not include the q0 point, whereas it does in chimat
       ! The main program assumes that nq and q(3,:) do not include q0
       ! Here we amend nq and q from chimat to make them consistent with the epsmat format
       nq = nq-nq0
       allocate (q (3,nq), stat=allocerr)
       do iq = 1,nq
          q(:,iq) = q_dummy(:,iq+nq0)
       enddo
    endif
    ! Setup FFT dimensions
    ! The size of the epsilon inverse /chi matrix in the binary file (nmtx)
    ! varies by q-point.
    ! We initialize the maximum dimension, maxdim, to 0 and then loop over the
    ! q-points in the binary file to find the maximum ntmx value, which is then
    ! set as maxdim.
    maxnfft = 0.d0
    maxdim = 0

    allocate (isortg (ng), stat=allocerr)
    allocate (isortgi (ng), stat=allocerr)

    do iq = 1,nq
       if (calctype .ne. 4) then
          ! We are not looking at chi: read epsmat file
          read(11) ngq, nmtx, (isortg(ig),isortgi(ig),ig=1,ng) ! read ngq, nmtx, isortq, isortqi
          read(11) ! read ekin
          read(11) ! read the current q point
          do i = 1,nmtx
             read(11) ! read epsinv
          enddo
       else
          ! We are looking at chi: read chimat
          read(09)
          read(09) nmtx
          do i = 1,nmtx
             read(09) ! read chi
          enddo
       endif

       ! Update maxdim
       if (nmtx .gt. maxdim) maxdim = nmtx

       ! Since nmtx is q-dependent so are nfft(:).
       ! The size of data_rrp is determined by the maximum values of nfft so here we
       ! set up the FFT grid in order to determinine the max nfft
       ! values and correctly allocate the size of data_rrp

       ! max (absolute) dimensions of G vectors in x, y and z directions
       kmax(1) = 0
       kmax(2) = 0
       kmax(3) = 0

!!!!!!LOG
       !write(85,*) "Data for q",iq,"of",nq
       !write(85,*) "nmtx value is",nmtx
       !write(85,*) "isortq is",(isortq(i),i=1,ng)
!!!!!!

       do i = 1,nmtx
          if (calctype .ne. 4) then
             ! isortg defined for epsmat file
             ig = isortg(i) ! return corresponding index from list of unsorted Gvecs
          else
             ! isortg not defined for chimat file
             ig = i
          endif
          if (abs(kx(ig)).gt.kmax(1)) then
             kmax(1) = abs(kx(ig))
          endif
          if (abs(ky(ig)).gt.kmax(2)) then
             kmax(2) = abs(ky(ig))
          endif
          if (abs(kz(ig)).gt.kmax(3)) then
             kmax(3) = abs(kz(ig))
          endif
       enddo
       ! Set max dimensions of FFT grid in x, y, and z.
       ! Note that Gvec dimensions run from -k to +k.
       nmax(1) = 2*kmax(1) + 1
       nmax(2) = 2*kmax(2) + 1
       nmax(3) = 2*kmax(3) + 1
       call setup_fft_sizes(nmax,nfft,scale)

       ! Set max nfft sizes
       do i = 1,size(nfft)
          if (nfft(i) .gt. maxnfft(i)) maxnfft(i) = nfft(i)
       enddo
    enddo

    ! Close appropriate binary file
    if (calctype .ne. 4) then
       close(11)
    else
       close(09)
    endif
    deallocate (isortg)
    deallocate (isortgi)

  end subroutine binaryfile_setup

  !------------------------------
  ! DEPRECATED -> REMOVE ME!
  ! Read inverse dielectric / chi matrix from
  ! appropriate BerkeleyGW binary file
  subroutine binaryfile_readmatrix(is_q0,initial_read,calctype)
    integer, intent(in) :: calctype
    logical, intent(inout) :: initial_read
    logical :: is_q0

    if (calctype .ne. 4) then
       ! Read eps(0)mat if we are not using chi
       if (is_q0) then
          ! We are looking at q0 so read from eps0mat

       else
          ! We are looking at a non-q0 point so read from epsmat

       endif

    else
       ! Otherwise read chimat
    endif
  end subroutine binaryfile_readmatrix

end module input_m

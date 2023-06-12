!=========================================================================
! Module input
! Contains routines relating to the input file for the program. 
! Adapted from Berkeley GW inread subroutine.
!=========================================================================

#include "f_defs.h"

module input_m

  use params_m
  use fftw_m
  implicit none

  public  :: read_input_wrrp
  public  :: read_input_wfcint
  public  :: read_input_autocorr
  public  :: read_q_points
  public  :: read_qs_weights
  public  :: read_wcoul0
  public  :: read_cube
  public  :: read_wrrp
  public  :: read_wrrp3d
  public  :: binaryfile_init
  public  :: binaryfile_setup
  public  :: binaryfile_readmatrix

contains

!------------------------------
! Read the input file 'wrrp.inp'
  subroutine read_input_wrrp(a,b,alat,blat,vol,use_q0,use_ss,use_wavg,use_slab,calctype,calcname,sc_size)
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

    ! Set default values
    ! KN: Currently values in wrrp.inp must be set in correct order
    !     This should be changed to a keyword system
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

  end subroutine read_input_wrrp

!------------------------------
! Read the input file 'wfcint.inp' for program wfcint.x
  subroutine read_input_wfcint(nat,atomtype,atompos,bader)
    integer :: ioerr,allocerr,iat
    integer, intent(out) :: nat
    character, allocatable, intent(out) :: atomtype*2(:)
    real(DP), allocatable,intent(out) :: atompos(:,:),bader(:)
    
    ! Read input file "wfcint.inp"
    open(unit=73,file='wfcint.inp',form='formatted',status='old')

    ! Read number of atoms    
    read(73,*,iostat=ioerr) nat
    ! Allocate variable dimensions
    allocate(atomtype(nat),stat=allocerr)
    allocate(atompos(3,nat),stat=allocerr)
    allocate(bader(nat),stat=allocerr)
    ! Read data for the atoms
    ! Atoms are listed in Cartesian coordinates (bohr)
    do iat = 1,nat
      read(73,*) atomtype(iat),atompos(1,iat),atompos(2,iat),atompos(3,iat),bader(iat)
    enddo
    close(73)
  end subroutine read_input_wfcint

!------------------------------
! Read the input file 'autocorr.inp' for program autocorr.x
  subroutine read_input_autocorr(zcoor,nextrachg,chgpos)
    integer :: ioerr, allocerr, i 
    integer, intent(out) :: nextrachg
    integer, intent(out) :: chgpos(3)
    real(DP), intent(out) ::zcoor

    ! Read input file "autocorr.inp"
    open(unit=73,file='autocorr.inp',form='formatted',status='old')
    ! Read specified z coordinate (in bohr) of slab surface
    read(73,*,iostat=ioerr) zcoor
    ! Read position of the charge (in units of FFT indices). This is the
    ! position of the perturbation in the 3D wrrp.3d.x calculation, and should
    ! be the same value as that found in the wrrp.3d.x input file
    read(73,*,iostat=ioerr) (chgpos(i),i=1,3)
    ! Read number of ADDITIONAL charges to be included in the cell. These will
    ! be distributed randomly in the same z plane specified by chgpos
    read(73,*,iostat=ioerr) nextrachg
    close(73)
  end subroutine read_input_autocorr

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
! Reads a cube file for the purpose of storing a QE-generated wavefunction
  subroutine read_cube(fname,ngrid,dela,atomtype,atompos,cubedata)
    character :: dummyc
    integer :: ioerr,allocerr,iat,i,j,k,nat
    real(DP) :: dummyr,origin
    character, intent(in) :: fname*32
    integer, intent(out) :: ngrid(3)
    real(DP), intent(out) :: dela(3,3) ! deltas of the lattice vectors, a()
    integer, allocatable, intent(out) :: atomtype(:)
    real(DP), allocatable, intent(out) :: cubedata(:),atompos(:,:)

    ! Open specified cube file
    open(66,file=fname,form='formatted',status='old',iostat=ioerr)
    
    ! Read data from cube file
    ! KN: Note that I'm assuming all units in bohr (i.e. ngrid values > 0)
    read(66,*) dummyc    ! ignore 1st header line
    read(66,*) dummyc    ! ignore 2nd header line
    read(66,*) nat, origin      ! number atoms and origin of cell
    read(66,*) ngrid(1),dela(:,1)       ! number of grid points and delta of a1
    read(66,*) ngrid(2),dela(:,2)       ! number of grid points and delta of a2
    read(66,*) ngrid(3),dela(:,3)       ! number of grid points and delta of a3
    ! read atom info into dummy variables
    allocate(atomtype (nat),stat=allocerr)
    allocate(atompos (3,nat),stat=allocerr)
    ! Read atom information (type and coordinates) for each atom
    do iat = 1,nat
      read(66,*) atomtype(iat),dummyr,atompos(:,iat)
    enddo
    
    allocate(cubedata (product(ngrid)), stat=allocerr)
    ! Read the cube data: outer loop corresponds to a1, middle to a2, and inner
    ! to a3
    read(66,*) cubedata
!    do i = 1,ngrid(1)
!      do j = 1,ngrid(2)
!        do k = 1,ngrid(3)
!          read(66,*) cubedata(
!        enddo
!      enddo
!    enddo
    close(66)
  end subroutine read_cube

!------------------------------
! Reads 6D wrrp.x output file
! KN: Note that this is presently a binary file
  subroutine read_wrrp(calctype,alat,a,nfft,sc_size,wrrpdata)
    integer :: i,j,ioerr,allocerr
    character, intent(out) :: calctype*12
    integer, intent(out) :: nfft(6), sc_size(3)
    real(DP), intent(out) :: alat,a(3,3)
    complex(DPC), allocatable, intent(out) :: wrrpdata(:,:,:,:,:,:)

    ! Open wrrp.x output file
    ! KN: This is hard-coded as "data_rrp.bin"
    open(89,file="data_rrp.bin",form='unformatted',status='old',iostat=ioerr)

    ! Read data from wrrp.x output file
    read(89)                            ! dummy string
    read(89) calctype                   ! wrrp.x calculation type
    read(89)                            ! dummy string
    do i = 1,3
      read(89) (a(i,j),j=1,3)           ! lattice vectors in alat units
    enddo
    read(89) alat                       ! alat in bohr
    read(89)                            ! dummy string
    read(89) (nfft(i),i=1,6)            ! FFT grid size
    read(89)                            ! dummy string
    read(89) (sc_size(i),i=1,3)         ! super cell size
    read(89)                            ! blank line

    allocate(wrrpdata(sc_size(1)*nfft(1),sc_size(2)*nfft(2),sc_size(3)*nfft(3),&
            sc_size(1)*nfft(4),sc_size(2)*nfft(5),sc_size(3)*nfft(6)),stat=allocerr)
    ! Read 6D data from wrrp.x output file    
    read(89) wrrpdata
    close(89)
  end subroutine read_wrrp

!------------------------------
! Reads 3D wrrp.x output file
! KN: Note that this is presently a binary file
  subroutine read_wrrp3d(calctype,alat,a,nfft,sc_size,wrrpdata)
    integer :: i,j,ioerr,allocerr
    character, intent(out) :: calctype*12
    integer, intent(out) :: nfft(3), sc_size(3)
    real(DP), intent(out) :: alat,a(3,3)
    complex(DPC), allocatable, intent(out) :: wrrpdata(:,:,:)

    ! Open wrrp.x output file
    ! KN: This is hard-coded as "data_rrp.bin"
    open(89,file="data_rrp.bin",form='unformatted',status='old',iostat=ioerr)

    ! Read data from wrrp.x output file
    read(89)                            ! dummy string
    read(89) calctype                   ! wrrp.x calculation type
    read(89)                            ! dummy string
    do i = 1,3
      read(89) (a(i,j),j=1,3)           ! lattice vectors in alat units
    enddo
    read(89) alat                       ! alat in bohr
    read(89)                            ! dummy string
    read(89) (nfft(i),i=1,3)            ! FFT grid size
    read(89)                            ! dummy string
    read(89) (sc_size(i),i=1,3)         ! super cell size
    read(89)                            ! blank line

    allocate(wrrpdata(sc_size(1)*nfft(1),sc_size(2)*nfft(2),sc_size(3)*nfft(3)),stat=allocerr)
    ! Read 6D data from wrrp.x output file    
    read(89) wrrpdata
    close(89)
  end subroutine read_wrrp3d

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

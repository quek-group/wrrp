!=========================================================================
! Calculates W(r,r') from BerkeleyGW static inverse dialectric matrix
! file. We first convert eps-1{GG'}(q;0) to real space eps-1(r,r') via FFTW,
! then calculate v(r,r') and perform the integral in real space.
! The reciprocal space eps-1 matrix is read from the epsmat binary file
! produced by the "Epsilon" BerkeleyGW code.
! Also capable of reading BerkeleyGW chimat files and printing chi(r,r').
!
! This code was written at the National University of Singapore (NUS)
! v0.1 K Noori
! v0.2 K Noori, N Cheng, FY Xuan
!=========================================================================

#include "f_defs.h"

program wrrp

  use params_m
  use fftw_m
  use input_m
  use supercell_m
  use v_m
  implicit none

  character :: dummy, calcname*12, date*8, time*10
  integer :: nq, nqs, ngq, nqfbz, nqstar, ng, nmtx, ig, igp, iq, iqf, iqs, nr, ir, irp, ifbz, i, j, ioerr, allocerr, maxdim, calctype

  ! XFY: variable descriptions
  ! nrq: number of q point in fz.
  ! ncouls: number of ig by eps cutoff.
  ! epstemp: eps matrix for q in fz.
  ! rqin: fz q point
  ! epsqpt: irbz q point for a corresponding fz q. indrq: fz -> irbz
  ! ind: g map; index of G_1 for a given G(ig). ph: phase factor exp(-i G tao)
  ! rq = r(q) + G(R) where rq in fz while q in irbz
  ! G_1 = R^(-1)( G+G(R) )
  integer :: irq, irq2, irq_, iscreen, nrq, ncouls, ncoulsm
  complex(DPC), allocatable :: ph(:,:)
  SCALAR, allocatable :: epstemp(:,:)
  integer, allocatable :: ind(:,:), indrq(:), isortg_fzeps(:,:)
  real(DP), allocatable :: rqin(:,:), epsqpt(:,:), q_fzeps(:,:)

  ! Nic modifications for supercell: 14.08.2018
  ! nr_sc     : total number of real space grid points
  ! sc_size   : User input for supercell size
  ! r_uc,rp_uc: real space grid points in the unit cell
  ! sc2uc_idx : mapping of fft index in the unit cell to the supercell
  ! scvec2real: mapping of supercell index to real space vector
  ! sc_a      : lattice vectors in the supercell
  ! sc_rcart  : real space grid points in the supercell
  ! vinv_gp   :
  ! sc_fftbox : fftbox for the supercell, with dimensions
  !             [sc_size(1)*maxnfft(1),...]
  integer :: k, l, m, n, nr_sc, sc_size(3)
  integer, allocatable :: r_uc(:), rp_uc(:), sc2uc_idx(:), scvec2real(:,:)
  real(DP) :: sc_a(3,3)
  real(DP), allocatable :: sc_rcart(:,:), vinv_gp(:,:)
  complex(DPC), allocatable :: sc_fftbox(:,:,:,:,:,:)
  SCALAR :: vcg, vcgp

  ! KN original variables
  integer :: x, y, z, xp, yp, zp
  integer :: kmax(3), nmax(3), nfft(6), maxnfft(6)
  logical :: use_q0, use_ss, use_wavg, use_slab, initial_read, use_symm
  real(DP) :: scale, a(3,3), b(3,3), alat, blat, vcell, r(3), rp(3), rtot(3), q0(3), qeff(3), qzero(3), qcart(3), rq, vcoul0, wcoul0r, wcoul0i, subcut, eking, ekingp
  complex(DPC) :: qfac, tmpsum
  SCALAR :: vcoul, wcoul0
  integer, allocatable :: isortg(:), isortgi(:), kx(:), ky(:), kz(:), glist(:,:), rvec2real(:,:), qmapi(:), qmapf(:)
  real(DP), allocatable :: q(:,:), qs(:,:), qfbz(:,:), qstar(:,:), wqs(:), ekin(:), rcart(:,:), v_gp(:,:)
  complex(DPC), allocatable :: fftbox(:,:,:,:,:,:), data_rrp(:,:,:,:,:,:), planavg(:)
  SCALAR, allocatable :: matrix_ggp(:,:), data_ggp(:,:), data_ggp1d(:)

  !==============================================================================
  ! Setup calculation
  write(6,*) repeat("=",79)
  call date_and_time(DATE=date,TIME=time)
  write(6,*) "wrrp.x -- A simple little program..."
  write(6,*) "version: beta"
  write(6,*) "start time: ",date," at ",time
  write(6,*) repeat("=",79)
#ifdef CPLX
  write(6,*) "We are running in COMPLEX mode."
#else
  write(6,*) "We are running in REAL mode."
#endif
  !x: open file which contains gmap ind(ig,irq) ph(ig,irq)
  open(19,file='gmapdata',form='unformatted',status='old')
  !x isortg contains all the isortg for each q in fz
  open(18,file='isortg',form='unformatted',status='old')
  ! Open log file
  open(unit=85,file="wrrp.log",form='formatted',status='replace',iostat=ioerr)

  ! Read input file, hardcoded as 'wrrp.inp'
  ! Nic: added new variable sc_size in wrrp.inp
  call read_input(a,b,alat,blat,vcell,use_q0,use_ss,use_wavg,use_slab,use_symm,calctype,calcname,sc_size)

  ! Convert cell volume from alat^3 to au^3 units
  vcell = alat**3 * vcell

  ! Output status of input read
  if (calctype .ne. 4) then
     write(6,*) "We will read a BerkeleyGW epsmat file."
     if (use_q0) then
        write(6,*) "We will read a BerkeleyGW eps0mat file."
     endif
  else
     write(6,*) "We will read a BerkeleyGW chimat file."
     if (use_q0) then
        write(6,*) "We will read a BerkeleyGW chi0mat file."
     endif
  endif

  ! If subsampling is enabled we read weights for subsampled q-pts
  if (use_ss) then
     write(6,*) "We are using the NNS scheme for q-->0."
     call read_qs_weights(nqs,wqs)
     ! Normalize weights
     wqs = wqs / sum(wqs)
     ! Determine cutoff for subsampled points
     ! Equivalent to Ekin of smallest G, |G|^2, in periodic direction
     ! In practice we ASSUME nonperiodic direction to be z and so subcut is
     ! min[Ekin(1,0,0),Ekin(0,1,0)]
     subcut = min(norm2(blat*b(:,1))**2,norm2(blat*b(:,2))**2)
     write(6,902) "Cutoff for subsampled q-points is",subcut,"Ry."
  endif

  ! If NOT using subsampling and if use_wavg is .true. in wrrp.inp then read
  ! vcoul0 and wcoul0 from wrrp.wcoul0.
  ! These values are generated by BerkeleyGW.
  ! If use_wavg .false. then set vcoul0 = wcoul0 = 0
  if (.not. use_ss) then
     if (use_wavg) then
        write(6,*) "Reading V and W head elements from wrrp.wcoul0..."
        call read_wcoul0(iscreen,vcoul0,wcoul0)
        !      open(65,file='wrrp.wcoul0',form='formatted',status='old',iostat=ioerr)
        !!      open(unit=121,file="slab.log",form='formatted',status='replace')
        !      read(65,*) dummy, iscreen ! BGW screening type
        !      read(65,*) !q0 vector
        !      read(65,*) !G0
        !      read(65,*) dummy, vcoul0
        !#ifdef CPLX
        !      read(65,*) dummy, wcoul0r, wcoul0i
        !      wcoul0 = cmplx(wcoul0r,wcoul0i)
        !#else
        !      read(65,*) dummy, wcoul0
        !#endif
        write(6,'(a9,f12.6)') "vcoul0: ",vcoul0
        write(6,'(a9,f12.6)') "wcoul0: ",wcoul0
        close(65)
     else
        write(6,*) "Setting V and W head elements to 0..."
        vcoul0 = 0.d0
        wcoul0 = 0.d0
        write(6,'(a9,f12.6)') "vcoul0: ",vcoul0
        write(6,'(a9,f12.6)') "wcoul0: ",wcoul0
     endif
  endif

  !==============================================================================
  ! Read q-point input files, hardcoded as 'wrrp.qibz' and 'wrrp.qfbz'
  ! KN: Filename is hardcoded and must be present for now
  write(6,*) "Reading q-point information from input..."
  call read_q_points(qfbz,nqfbz,qmapf,qmapi)

  !==============================================================================
  ! Read binary file and setup calculation
  write(6,*) "Reading binary file(s) and setting up calculation..."

  ! Initial file read to extract ng and nq
  call binaryfile_setup(nq,nqs,ng,kx,ky,kz,q,qs,q0,use_q0,use_ss,nfft,maxnfft,calctype)

!!!!!LOG
  write(85,900) "The lattice constant alat (in bohr) is",alat
  write(85,*) "The direct lattice vectors (in alat) are:"
  write(85,904) ((a(i,j),i=1,3),j=1,3)
  write(85,900) "The cell volume in au^3 is",vcell
  write(85,*) repeat("=",79)
  write(85,900) "The reciprocal lattice constant blat (in 2*pi/alat) is",blat
  write(85,*) "The reciprocal lattice vectors (in blat) are:"
  write(85,904) ((b(i,j),i=1,3),j=1,3)
  write(85,*) repeat("=",79)
  write(85,*) "There are",nqfbz,"q-points in the full BZ."
  if (use_q0) then
     if (use_ss) then
        write(85,*) "There are",nqs,"subsampled q-points (qs) for q-->0:"
        write(85,904) ((qs(i,iq),i=1,3),iq=1,nqs)
        write(85,902) "The dielectric matrix cutoff for qs-points is",subcut,"Ry."
     else
        write(85,*) "There is one q0 point:"
        write(85,903) q0
     endif
  endif
  write(85,*) "There are",nq,"non-q0 q-points in the irreducible BZ:"
  write(85,904) ((q(i,iq),i=1,3),iq=1,nq)
  write(85,*) repeat("=",79)
  write(85,*) "There are",ng,"G vectors."
  !write(85,*) "The G vectors are:"
  !write(85,*) (kx(i),ky(i),kz(i),i=1,ng)
  write(85,*) repeat("=",79)
900 format(1x,a,f16.9)
901 format(1x,a,3f16.9)
902 format(1x,a,f16.9,1x,a)
903 format(f16.9)
904 format(3f16.9)
!!!!!LOG

  ! Store unsorted G vectors in a master list
  allocate (glist (3,ng), STAT=allocerr)
  do i = 1, ng
     glist(1,i) = kx(i)
     glist(2,i) = ky(i)
     glist(3,i) = kz(i)
  enddo

  !TEST: manually change maxnfft
  !maxnfft(1) = 16
  !maxnfft(2) = 16
  !!maxnfft(3) = 24
  !maxnfft(4) = 16
  !maxnfft(5) = 16
  !maxnfft(6) = 24
!!!!!

  scale = 1.0d0/product(maxnfft(1:6))
  !write(6,'(a,e14.6)') "The FFT scaling factor is ",scale
  if (use_q0) then
     if (use_ss) then
        write(6,'(1x,a,i6,1x,a)') "There are",nqs,"subsampled q-points for q-->0."
     else
        write(6,*) "There is a q0 point."
     endif
  endif
  write(6,'(1x,a,i6,1x,a)') "There are",nq,"non-q0 q-points in the irreducible BZ."
  write(6,'(1x,a,i6,1x,a)') "There are",nqfbz,"q-points in the full BZ."
  write(6,'(1x,a,f9.3,1x,a)') "The cell volume is",vcell,"au^3."
  write(6,*) "The maximum 6D FFT box size is:"
  write(6,'(i3,i3,i3,i3,i3,i3)') (maxnfft(i),i=1,6)

!!!!!LOG
  write(85,*) "The maximum FFT box size is:"
  write(85,'(i3,i3,i3,i3,i3,i3)') (maxnfft(i),i=1,6)
!!!!!LOG

  !  allocate (data_rrp (maxnfft(1),maxnfft(2),maxnfft(3),maxnfft(4),maxnfft(5),maxnfft(6)), STAT=allocerr)
  allocate (fftbox (maxnfft(1),maxnfft(2),maxnfft(3),maxnfft(4),maxnfft(5),maxnfft(6)), STAT=allocerr)

  ! Determine real space grid
  write(6,*) "Determining real-space grid..."
  allocate (rcart (3,maxnfft(1)*maxnfft(2)*maxnfft(3)))
  allocate (rvec2real (3,maxnfft(1)*maxnfft(2)*maxnfft(3)))
  nr = size(rcart,2)

  ! Nic modifications
  ! Allocation of different arrays needed to expand to supercell
  ! HELP: how can we improve memory usage? Memory usage is horrible
  nr_sc = nr*sc_size(1)*sc_size(2)*sc_size(3)
  allocate (sc2uc_idx(nr_sc))
  allocate (scvec2real(3,nr_sc))
  allocate (sc_rcart(3,nr_sc))
  allocate (r_uc(nr_sc))
  allocate (rp_uc(nr_sc))
  allocate (data_rrp (sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3),&
       sc_size(1)*maxnfft(4),sc_size(2)*maxnfft(5),sc_size(3)*maxnfft(6)), STAT=allocerr)
  allocate (sc_fftbox (sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3),&
       sc_size(1)*maxnfft(4),sc_size(2)*maxnfft(5),sc_size(3)*maxnfft(6)), STAT=allocerr)
  ! supercell.f90 used in get_real_grid
  call get_real_grid(alat,a,maxnfft,rcart,rvec2real,sc_size,sc2uc_idx, &
       scvec2real,sc_rcart,sc_a)
  !  call get_real_grid(alat, a, maxnfft, rcart, rvec2real)

  !==============================================================================
  ! W loop for q0 (i.e. q0 point)
  write(6,*) repeat("=",79)

  data_rrp = 0.d0
  qzero = 0.d0 ! Define a true q = 0 vector

  if (use_q0) then
     ! If use_q0 = .true. then read eps0mat/chi0mat file first
     write(6,*) "Entering ",calcname,"loop for q-->0."

     if (calctype .ne. 4) then
        ! We are not calculating chi; read eps0mat file
        open(12,file='eps0mat',form='unformatted',status='old',iostat=ioerr)

        if (use_ss) then
           ! We are using NNS subsampling for q-->0:
           ! eps0mat contains nqs q-points, qs

           ! Read the header info
           read(12)
           read(12)
           read(12)
           read(12)
           read(12)
           read(12)
           read(12)
           !        read(12) nqs
           !        backspace(12)
           !        allocate (qs (3,nqs),stat=allocerr)
           read(12) !nqs,((qs(i,iq),i=1,3),iq=1,nqs)
           read(12) !ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)

           ! Loop over qs points: note that q0 == qs(1)
           do iq = 1,nqs
              ! Read binary file for current qs-point. We allocate here since the
              ! matrix size, nmtx, is q-dependent.
              read(12) ngq
              backspace(12)
              allocate (isortg (ngq), stat=allocerr)
              allocate (isortgi (ngq), stat=allocerr)
              read(12) ngq, nmtx, (isortg(ig),isortgi(ig),ig=1,ngq)
              !write(6,*) "iq is",iq,"nmtx is",nmtx,"ngq is",ngq
              allocate (ekin (ngq), stat=allocerr)
              !          allocate (v_gp (nmtx,nmtx), stat=allocerr)
              allocate (matrix_ggp (nmtx,nmtx), stat=allocerr)
              ! Nic
              !          allocate (vinv_gp (nmtx,nmtx), stat=allocerr)
              !write(6,*) "size of eps is",size(matrix_ggp)
              ! The size of data_ggp will be set by the rank of epsinv(q0) (where q0
              ! = qs(1)) since the remaining qs-points only contribute to the neck
              ! elements (i.e. small G-vectors). We therefore only allocate and
              ! initialize these values for q0.
              if (iq .eq. 1) then
                 allocate (data_ggp (nmtx,nmtx), stat=allocerr)
                 allocate (data_ggp1d (nmtx*nmtx), stat=allocerr)
                 data_ggp = 0.d0
                 data_ggp1d = 0.d0
              endif
              !          v_gp = 0.d0

              read(12) (ekin(ig),ig=1,ngq)
              read(12) ! read the current q point
              do igp = 1,nmtx
                 ! read the epsinv matrix column-by-column; G vectors sorted by KE
                 read(12) (matrix_ggp(ig,igp),ig=1,nmtx)
              enddo
              !write(6,*) "I've read epsmat and nmtx is",nmtx
              !          ! Construct asymmetric V matrix
              !          do igp = 1,nmtx ! G' loop
              !            !call vcoul_g(slab,glist(:,isortg(igp)),q(:,iq),b,blat,vcoul0,vcoul)
              !            call vcoul_g(use_slab,glist(:,isortg(igp)),qs(:,iq),b,blat,vcoul0,vcoul)
              !            do ig = 1,nmtx ! G loop
              !              if (ig .eq. igp) then
              !                v_gp(ig,igp) = vcoul ! diagonal elements (G=G')
              !              else
              !                v_gp(ig,igp) = 0.d0  ! off-diagonal elements are 0
              !              endif
              !            enddo
              !          enddo

              ! Calculate specified property: note that only W or Vscr
              ! should be computed when using subsampling.
              select case (calctype)
              case (0)
                 write(6,*) "Calculating W(G,G') for qs-point",iq,"of",nqs,"with a weight of ",wqs(iq)
              case (1)
                 write(6,*) "Calculating Vscr(G,G') for qs-point",iq,"of",nqs,"with a weight of ",wqs(iq)
              case (5)
                 write(6,*) "Calculating rhoind(G,G') for qs-point",iq,"of",nqs,"with a weight of ",wqs(iq)
              end select
              write(6,'(3f10.6)') qs(:,iq)

              ! Set effective q that is passed for calculation of V
              qeff = qs(:,iq) ! use the current qs-point

              ! Nic according to xuan's suggestions:
              ! Compute interacting chi and save to matrix_ggp first
              !          if (calctype .eq. 5) then
              !            ! Construct asymmetric V matrix
              !            do igp = 1,nmtx ! G' loop
              !              !call vcoul_g(slab,glist(:,isortg(igp)),q(:,iq),b,blat,vcoul0,vcoul)
              !              call vcoul_g(use_slab,glist(:,isortg(igp)),qs(:,iq),b,blat,vcoul0,vcoul)
              !              do ig = 1,nmtx ! G loop
              !                if (ig .eq. igp) then
              !                  vinv_gp(ig,igp) = 1.d0/vcoul ! diagonal elements (G=G')
              !                else
              !                  vinv_gp(ig,igp) = 0.d0  ! off-diagonal elements are 0
              !                endif
              !              enddo
              !            enddo
              !
              !            ! matmul as usual
              !            matrix_ggp = matmul(vinv_gp,matrix_ggp) - vinv_gp !chi
              !          endif

              ! If |G|^2 AND |G'|^2 > subcut then we set q0=qs(1)
              ! and compute the property as normal.
              ! Otherwise we have a neck element and therefore perform the
              ! sum over qs. See eqs. 10-11 PRB 95, 035109 (2017)
              do igp = 1,nmtx ! G' loop

                 !            ekingp = matmul(glist(:,isortg(igp)),b)
                 ekingp = norm2(blat*matmul(glist(:,isortg(igp)),b))**2 ! KE of the current G'
                 ! compute asymmetric V(q+G')
                 call vcoul_g(use_slab,glist(:,isortg(igp)),qeff,b,blat,vcoul0,vcoul)
                 vcgp = vcoul
                 do ig = 1,nmtx ! G loop
                    call vcoul_g(use_slab,glist(:,isortg(igp)),qeff,b,blat,vcoul0,vcoul)
                    vcg = vcoul
                    !              ! fill asymmetric V matrix
                    !              if (ig .eq. igp) then
                    !                ! diagonal elements only
                    !                v_gp(ig,igp) = vcoul
                    !              endif
                    ! Check kinetic energy of the G and G' vectors
                    eking = norm2(blat*matmul(glist(:,isortg(ig)),b))**2 ! KE of the current G
                    if ((eking .ge. subcut) .and. (ekingp .ge. subcut)) then
                       ! non-neck element: only deal with q0 = qs(1)
                       if (iq .eq. 1) then
                          ! compute W
                          if (calctype .ne. 5) then
                             data_ggp(ig,igp) = matrix_ggp(ig,igp)*vcgp
                             ! Compute Vscr by substracting V when G = G'
                             if ((calctype .eq. 1) .and. (ig .eq. igp)) then
                                data_ggp(ig,igp) = data_ggp(ig,igp) - vcgp
                             endif
                          else
                             if (ig .eq. igp) then
                                data_ggp(ig,igp) = (matrix_ggp(ig,igp) - 1.0d0) / vcg
                                data_ggp(ig,igp) = data_ggp(ig,igp) * vcgp
                             else
                                data_ggp(ig,igp) = matrix_ggp(ig,igp)  / vcg
                                data_ggp(ig,igp) = data_ggp(ig,igp) * vcgp
                             endif
                          endif
                       endif
                    else
                       if (calctype .ne. 5) then
                          ! neck element:
                          ! compute weighted W and add to previous value (i.e. from other qs)
                          data_ggp(ig,igp) = data_ggp(ig,igp) + (wqs(iq)*matrix_ggp(ig,igp)*vcgp)
                          ! Compute Vscr by substracting weighted V when G = G'
                          if ((calctype .eq. 1) .and. (ig .eq. igp)) then
                             data_ggp(ig,igp) = data_ggp(ig,igp) - (wqs(iq)*vcgp)
                          endif
                       else
                          if (ig .eq. igp) then
                             data_ggp(ig,igp) = data_ggp(ig,igp) +&
                                  (wqs(iq)*(((matrix_ggp(ig,igp) - 1.0d0) / vcg) * vcgp ))
                          else
                             data_ggp(ig,igp) = data_ggp(ig,igp) +&
                                  (wqs(iq)*((matrix_ggp(ig,igp) / vcg) * vcgp ))
                          endif
                       endif
                    endif
                 enddo ! end ig loop
              enddo ! end igp loop

           enddo ! end iq loop

        else
           ! We are NOT using NNS subsampling for q-->0:
           ! eps0mat will only have one q-point, q0.
           ! We calcuate W(G,G';q=0) = eps^-1(G,G';q=q0)*v(G';q=0), replacing the
           ! head elements of W and V with wcoul0 and vcoul0, respectively.
           write(6,*) "Shifted q-point, q0, is:"
           write(6,'(3f10.6)') q0
           ! Read header information
           read(12)
           read(12)
           read(12)
           read(12)
           read(12)
           read(12)
           read(12)
           read(12) !nq,((q(i,iq),i=1,3),iq=1,nq)
           read(12) !ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)

           ! Read binary file for current q point
           read(12) ngq
           backspace(12)
           allocate (isortg (ngq), stat=allocerr)
           allocate (isortgi (ngq), stat=allocerr)
           read(12) ngq, nmtx, (isortg(ig),isortgi(ig),ig=1,ngq) ! nmtx is q-dependent
           allocate (ekin (ngq), stat=allocerr)
           allocate (v_gp (nmtx,nmtx), stat=allocerr)
           ! Nic
           allocate (vinv_gp (nmtx,nmtx), stat=allocerr)
           allocate (matrix_ggp (nmtx,nmtx), stat=allocerr)
           allocate (data_ggp (nmtx,nmtx), stat=allocerr)
           allocate (data_ggp1d (nmtx*nmtx), stat=allocerr)
           data_ggp = 0.d0
           data_ggp1d = 0.d0
           v_gp = 0.d0
           vinv_gp = 0.d0
           read(12) (ekin(ig),ig=1,ngq)
           read(12) ! read the current q point
           do igp = 1,nmtx
              ! read the epsinv matrix column-by-column; G vectors sorted by KE
              read(12) (matrix_ggp(ig,igp),ig=1,nmtx)
           enddo

           ! Set effective q that is passed for calculation of V
           if ((iscreen .eq. 2) .and. use_slab) then
              ! We have a (2D) truncated metal so we pass the q0 vector. Note that
              ! this is because, in such a case, epsilon^-1 and V are independent
              ! of q as q-->0 (see BGW paper and D Strubbe thesis).
              qeff = q0
           else
              ! Otherwise we directly pass the q = (0,0,0) vector
              !          qeff = qzero
              qeff = q0
           endif

           ! Construct asymmetric V matrix
           do igp = 1,nmtx ! G' loop
              !call vcoul_g(slab,glist(:,isortg(igp)),q(:,iq),b,blat,vcoul0,vcoul)
              !call vcoul_g(use_slab,glist(:,isortg(igp)),q0,b,blat,vcoul0,vcoul)
              call vcoul_g(use_slab,glist(:,isortg(igp)),qeff,b,blat,vcoul0,vcoul)
              do ig = 1,nmtx ! G loop
                 if (ig .eq. igp) then
                    ! G = G' : diagonal element = vcoul
                    v_gp(ig,igp) = vcoul
                    vinv_gp(ig,igp) = 1.d0/vcoul
                    ! Nic: problem if vcoul happen to be zero, replace with 1e-12 for
                    ! now
                    ! Can't be 1e+12... try vcoul0 for now
                    !              if (vcoul .eq. 0.d0) then
                    !                vinv_gp(ig,igp) = 0.d0
                    !              endif
                 else
                    ! G != G' : off-diagonal element = 0
                    v_gp(ig,igp) = 0.d0
                    vinv_gp(ig,igp) = 0.d0
                 endif
              enddo
           enddo

        endif

     else
        ! We are calculating chi; read chi0mat file
        open(10, file='chi0mat',form='unformatted',status='old',iostat=ioerr)
        ! Read header information
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)
        read(10)

        ! Read binary file for current q point
        read(10)
        read(10) nmtx
        allocate (matrix_ggp (nmtx,nmtx), stat=allocerr)
        allocate (data_ggp (nmtx,nmtx), stat=allocerr)
        allocate (data_ggp1d (nmtx*nmtx), stat=allocerr)
        do igp = 1,nmtx
           ! read the chi0mat matrix column-by-column
           read(10) (matrix_ggp(ig,igp),ig=1,nmtx)
        enddo
     endif

     ! Calculate specified property if not using NNS scheme.
     ! Calculation of NNS properties is done in iq loop above.
     if (.not. use_ss) then
        select case (calctype)
        case (0)
           write(6,*) "Calculating W(G,G')..."
           ! calculate W(G,G')
           data_ggp = matmul(matrix_ggp,v_gp)
           ! replace head element W(0,0; q=q0) with wcoul0
           write(6,*) "Replacing head element of W with wcoul0..."
           data_ggp(1,1) = wcoul0
        case(1)
           ! KN: FIX THIS: IS THIS THE CORRECT WAY TO REPLACE wcoul0 and vcoul0
           ! for Vscr? Should vcoul0 be included in v_gp during multiplication of
           ! eps^-1 and V?
           write(6,*) "Calculating Vscr(G,G')..."
           ! calculate W(G,G')
           data_ggp = matmul(matrix_ggp,v_gp)
           ! replace head element W(0,0; q=q0) with wcoul0
           write(6,*) "Replacing head element of W with wcoul0..."
           data_ggp(1,1) = wcoul0
           !v_gp(1,1) = vcoul0
           ! calculate Vscr(G,G') = W(G,G') - V(G')
           data_ggp = data_ggp - v_gp
        case(2)
           write(6,*) "Calculating V(G')..."
           ! return bare Coulomb potential, V
           data_ggp = v_gp
        case(3)
           write(6,*) "Returning epsinv(G,G')..."
           ! return epsilon inverse
           data_ggp = matrix_ggp
        case(4)
           write(6,*) "Returning chi(G,G')..."
           ! return chi
           data_ggp = matrix_ggp
        case(5)
           write(6,*) "Calculating rhoind(G,G')..."
           ! calculate rhoind(G,G')
           ! First calculate interacting chi (chi)
           ! Eq. (11) PRB 35 5585 (1987)
           ! epsinv(g,gp) = 1 + v(gp)*chi(g,gp)
           ! chi(g,gp) = v^-1(gp) * (epsinv(g,gp) - 1)
           data_ggp = matmul(vinv_gp,matrix_ggp) - vinv_gp !chi
           ! Eq. (1) PRB 35 5585 (1987)
           ! Taking the external perturbation to be a (positive)
           ! test charge,
           ! rhoind(g,gp) = chi(g,gp)*v(gp)
           data_ggp = matmul(data_ggp,v_gp) !rhoind
        end select
     endif

     ! Put W_GG(q) data into FFT box in preparation for FFTW
     ! Store W_GG(q) matrix in 1D array suitable for FFTW.
     ! Stored as: [G1G'1, G1G'2, ..., G1G'N, G2G'1, G2G'2, ..., G2G'N, ...,
     ! GNG'1, GNG'2, ..., GNG'N]
     i = 1
     do ig = 1,nmtx
        do igp = 1,nmtx
           data_ggp1d(i) = data_ggp(ig,igp)
           i = i + 1
        enddo
     enddo

     ! Put data into box assuming G -> G and G' -> -G'
     if (calctype .ne. 4) then
        ! We are reading eps0mat and need isortg
        call put_into_fftbox(data_ggp1d, nmtx, ng, glist, fftbox, maxnfft, isortg)
     else
        ! We are reading chi0mat and don't need isortg
        call put_into_fftbox(data_ggp1d, nmtx, ng, glist, fftbox, maxnfft)
     endif

     !write(6,*) fftbox
     ! Do the backward FFT
     call do_fft(fftbox,maxnfft,1) ! -1 for forward FFT and 1 for backward FFT

     ! Determine scaling factor for current q point
     ! From the Hybertsen/Louie paper we note that the inverse FFT of W_GG'(q)
     ! must be scaled by the factor exp[iq.(r-r')] for every q point.
     write(6,*) "Scaling ",calcname,"for q-->0..."
     ! Convert q to Cartesian 2pi/a units
     !qcart = 0.d0

     ! q is 0 for the purpose of q-pt scaling
     ! KN: even for truncated metals, right?
     qcart = qzero
     ! Nic: maybe qzero is why it becomes periodic?
     !    qcart = 0.d0
     !    qcart = q0(1)*b(:,1) + q0(2)*b(:,2) + q0(3)*b(:,3)

     ! Convert q to Cartesian 1/au units
     !    qcart = qcart * blat

     !write(6,*) fftbox(:,1,1,1,1,1)
     call get_sc_fftbox(fftbox,nr_sc,scvec2real,rvec2real,sc2uc_idx,sc_fftbox)
     !write(6,*) fftbox(:,1,1,1,1,1)
     !write(6,*) sc_fftbox(:,1,1,1,1,1)
     write(6,*) "Copy to fftbox success"
     ! Loop over r and r' in order to multiply by the required exp[q(r-r')]
     ! factor
     !$omp parallel do default(private) shared(sc_rcart,scvec2real,qcart,sc_fftbox,nr_sc)
     do ir = 1,nr_sc ! r loop
        do irp = 1,nr_sc ! r' loop
           qfac = (0.d0,0.d0)
           r = sc_rcart(:,ir) ! r
           rp = sc_rcart(:,irp) ! r'
           rtot = r - rp
           rq = dot_product(qcart,rtot)
           qfac = exp(I_D*rq)
           ! Now scale W(r,r') by the corresponding factor for the current q
           sc_fftbox(scvec2real(1,ir),scvec2real(2,ir),scvec2real(3,ir),scvec2real(1,irp),scvec2real(2,irp),scvec2real(3,irp)) = &
                sc_fftbox(scvec2real(1,ir),scvec2real(2,ir),scvec2real(3,ir),scvec2real(1,irp),scvec2real(2,irp),scvec2real(3,irp))  * &
                qfac
        enddo
     enddo
     !$omp end parallel do
     ! Add W_rr'(q0) to previous W_rr'
     ! Nic: fftbox -> sc_fftbox
     data_rrp = data_rrp + sc_fftbox

     ! Clean up
     if (calctype .ne. 4) then
        close(12)
        deallocate (isortg)
        deallocate (isortgi)
        deallocate (ekin)
        if (.not. use_ss) then
           deallocate (v_gp) ! not used in subsampling mode
           ! Nic
           deallocate (vinv_gp)
        endif
     else
        close(10)
     endif
     deallocate (matrix_ggp)
     deallocate (data_ggp)
     deallocate (data_ggp1d)

     write(6,*) repeat("=",79)
  endif

  !data_rrp = 0.d0

  !==============================================================================
  ! Nic: If we are not using symmetry, skip this 
  if (use_symm) then
     ! XFY Modifications for enabling symmetry
     !x: read in ind(ig,irq) ph(ig,irq) generates by gmap using BGW sigma with one specific point only.
     !x number of q points in fz
     read(19) nrq
     !x number of g vectors due to wave function cutoff
     read(19) ngq
     allocate (isortg_fzeps(ngq,nrq))
     allocate (q_fzeps(3,nrq))

     allocate (rqin(3,nrq))
     allocate (epsqpt(3,nrq))
     write(6,*)
     write(6,*) "xuan: nrq = ", nrq

     irq = 1
     read(19) ncouls
     ncoulsm = ncouls + 200
     allocate (ind(ncoulsm,nrq))
     allocate (ph(ncoulsm,nrq))
     allocate (indrq(nrq))

     irq = 1
     !x G -> G_1 for each q in fz
     read(19) (ind(ig,irq),ig=1,ncouls)
     read(19) (ph(ig,irq),ig=1,ncouls)
     !x q points in fz together with gmap
     read(19) (rqin(i,irq),i=1,3)
     read(19) (indrq(irq))
     !x q points in irbz corresponding to each q in fz
     read(19) (epsqpt(i,irq),i=1,3)

     do irq = 2,nrq
        !x isortg readin from file for each q in fz
        read(18) (isortg_fzeps(ig,irq),ig=1,ngq)
        !x fz q points together with isortg
        read(18) (q_fzeps(i,irq),i=1,3)

        read(19) ncouls
        read(19) (ind(ig,irq),ig=1,ncouls)
        read(19) (ph(ig,irq),ig=1,ncouls)
        read(19) (rqin(i,irq),i=1,3)

        read(19) (indrq(irq))
        read(19) (epsqpt(i,irq),i=1,3)
     enddo
  endif

  !==============================================================================
  ! Main loop for all other q points

  initial_read = .true.

  do iq = 1,nq ! irreducible q loop
     write(6,'(a9,a13,a36,i4,a3,i4,a1)') "Entering",calcname,"loop for non-q0 irreducible q-point",iq,"of",nq,":"
     write(6,'(3f10.6)') (q(j,iq),j=1,3)

     if (calctype .ne.4) then
        ! We are not calculating chi; read epsmat file

        ! Read header info
        if (initial_read) then
           open(11,file='epsmat',form='unformatted',status='old',iostat=ioerr)
           read(11)
           read(11)
           read(11)
           read(11)
           read(11)
           read(11)
           read(11)
           read(11) !nq,((q(i,iq),i=1,3),iq=1,nq)
           read(11) !ng,(kx(ig),ky(ig),kz(ig),ig=1,ng)
           initial_read = .false.
        endif
        ! Read q-dependent info
        read(11) ngq
        backspace(11)
        allocate (isortg (ngq), stat=allocerr)
        allocate (isortgi (ngq), stat=allocerr)
        read(11) ngq, nmtx, (isortg(ig),isortgi(ig),ig=1,ngq) ! nmtx is q-dependent
        allocate (ekin (ngq), stat=allocerr)
        allocate (matrix_ggp (nmtx,nmtx), stat=allocerr)
        !x
        allocate (epstemp (nmtx,nmtx), stat=allocerr)
        !Nic
        allocate (vinv_gp (nmtx,nmtx), stat=allocerr)
        allocate (v_gp (nmtx,nmtx), stat=allocerr)
        allocate (data_ggp (nmtx,nmtx), stat=allocerr)
        allocate (data_ggp1d (nmtx*nmtx), stat=allocerr)
        read(11) (ekin(ig),ig=1,ngq)
        read(11) ! read the current q point
        ! read the epsinv matrix column-by-column; G vectors sorted by KE
        do igp = 1,nmtx
           read(11) (matrix_ggp(ig,igp),ig=1,nmtx)
        enddo

        ! Nic: copy to epstemp
        epstemp = matrix_ggp

        write(6,*) "Unfolding current q-point..."
        ! If q0 then qmapi has length (nq+1), with qmapi(1) = qmapi_q0, so we must add 1 to qmapi index
        if (use_q0) then
           ifbz = qmapi(iq+1)
        else
           ifbz = qmapi(iq)
        endif
        ! Find unfolded q-points corresponding to current q-point in iBZ
        nqstar = 1
        do iqf = 1,nqfbz
           if (qmapf(iqf) .eq. ifbz) then
              nqstar = nqstar + 1
           endif
        enddo
        write(6,'(a10,i4,a32)') "There are",nqstar,"q-points in the current q star:"
        allocate(qstar (3,nqstar), stat=allocerr)
        ! Set the first q-point in the unfolded star to be the current folded
        ! q-point
        qstar(:,1) = q(:,iq)
        ! Add the remaning unfolded q-points, if present, to the star
        if (nqstar .gt. 1) then
           i = 1
           do iqf = 1,nqfbz
              if (qmapf(iqf) .eq. ifbz) then
                 i = i + 1
                 qstar(:,i) = qfbz(:,iqf)
              endif
           enddo
        endif
        do iqs = 1,nqstar
           write(6,'(3f10.6)') qstar(:,iqs)
        enddo

        ! Loop over the q-points in the current star
        do iqs = 1,nqstar
           data_ggp = 0.d0
           data_ggp1d = 0.d0
           v_gp = 0.d0
           vinv_gp = 0.d0
           ! Nic: skip if symmetry not used
           !x: using BerkeleyGW sigma code (one point only, since subgroup sysmmetry operations of Gamma poin
           !   give the smallest irbz) to generate ind(ig,irq) ph(ig,irq ) which maps iqs --> iq
           !   search which q in fz this iqs corresponds to, then save it to irq
           if (use_symm) then
              irq_loop: do irq_ = 1,nrq
                 !x  find which index of gmap produce the current q point in qstar loop
                 if ( all( abs( rqin(1:3,irq_) - qstar(1:3,iqs)) .lt. TOL_Small )) then
                    !                write(6,*) "1st if in"
                    if (all( abs( epsqpt(1:3,irq_) - q(1:3,iq)) .lt. TOL_Small)) then
                       !                  write(6,*) "2nd if in"
                       irq=irq_
                       exit irq_loop
                    endif
                 endif
              enddo irq_loop ! irq_

              irq2_loop: do irq_ = 2,nrq
                 !x  find which index of isortg produce the current q point in qstar loop
                 if ( all( abs( q_fzeps(1:3,irq_) - qstar(1:3,iqs)) .lt. TOL_Small )) then
                    irq2=irq_
                    exit irq2_loop
                 endif
              enddo irq2_loop ! irq_
              !x use isortg readin from file
              do igp = 1,ngq
                 isortg(igp)=isortg_fzeps(igp,irq2)
              enddo
           endif
              !x Construct epsmat at iqs (fz) using iq (irz)

           ! Construct asymmetric v(q+G') matrix
           do igp = 1,nmtx ! G' loop
              call vcoul_g(use_slab,glist(:,isortg(igp)),qstar(:,iqs),b,blat,vcoul0,vcoul)
              do ig = 1,nmtx ! G loop
                 if (ig .eq. igp) then
                    v_gp(ig,igp) = vcoul ! diagonal elements (G=G')
                    vinv_gp(ig,igp) = 1.0d0/vcoul
                 else
                    v_gp(ig,igp) = 0.d0  ! off-diagonal elements are 0
                    vinv_gp(ig,igp) = 0.d0
                 endif

                 if (use_symm) then
                    !x            q map for epsilon matrix
                    epstemp(ig,igp)=ph(ig,irq)*conjg(ph(igp,irq))*matrix_ggp(ind(ig,irq),ind(igp,irq))
                 endif 
              enddo ! end G loop
           enddo ! end G' loop

           ! Calculate specified property
           if (calctype .eq. 0) then
              write(6,'(a35,i4,a3,i4,a4)') "Calculating W(G,G') for unfolded q",iqs,"of",nqstar,"..."
              ! calculate W(G,G')
              data_ggp = matmul(epstemp,v_gp)
           elseif (calctype .eq. 1) then
              write(6,'(a38,i4,a3,i4,a4)') "Calculating Vscr(G,G') for unfolded q",iqs,"of",nqstar,"..."
              ! calculate Vscr(G,G') = W(G,G') - V(G')
              data_ggp = matmul(epstemp,v_gp) - v_gp
           elseif (calctype .eq. 2) then
              write(6,'(a33,i4,a3,i4,a4)') "Calculating V(G') for unfolded q",iqs,"of",nqstar,"..."
              ! return bare Coulomb potential, V
              data_ggp = v_gp
           elseif (calctype .eq. 3) then
              write(6,'(a40,i4,a3,i4,a4)') "Calculating epsinv(G,G') for unfolded q",iqs,"of",nqstar,"..."
              ! return epsilon inverse
              data_ggp = epstemp
              ! Nic
           elseif (calctype .eq. 5) then
              write(6,'(a40,i4,a3,i4,a4)') "Calculating fullchi(G,G') for unfolded q",iqs,"of",nqstar,"..."
              ! First calculate interacting chi (chi)
              ! Eq. (11) PRB 35 5585 (1987)
              ! epsinv(g,gp) = 1 + v(gp)*chi(g,gp)
              ! chi(g,gp) = v^-1(gp) * (epsinv(g,gp) - 1)
              data_ggp = matmul(vinv_gp,epstemp) - vinv_gp !fullchi
              ! Eq. (1) PRB 35 5585 (1987)
              ! Taking the external perturbation to be a (positive)
              ! test charge,
              ! rhoind(g,gp) = chi(g,gp)*v(gp)
              data_ggp = matmul(data_ggp,v_gp)
           else
              ! There is an error
           endif

           ! Put W_GG(q) data into FFT box in preparation for FFTW
           ! Store W_GG(q) matrix in 1D array suitable for FFTW.
           ! Stored as: [G1G'1, G1G'2, ..., G1G'N, G2G'1, G2G'2, ..., G2G'N, ...,
           ! GNG'1, GNG'2, ..., GNG'N]
           i = 1
           do ig = 1,nmtx
              do igp = 1,nmtx
                 !w_g1d(i) = w_g(ig,igp)
                 data_ggp1d(i) = data_ggp(ig,igp)
                 i = i + 1
              enddo
           enddo

           ! Put data into FFT box
           ! Note that this routine places the data assuming G -> G and G' -> -G'
           ! since we have W_rr'(q) = exp[i(q+G)r].W_GG'(q).exp[-i(q+G')r']
           call put_into_fftbox(data_ggp1d, nmtx, ng, glist, fftbox, maxnfft, isortg)

           ! Do the backward FFT
           call do_fft(fftbox,maxnfft,1) ! -1 for forward FFT and 1 for backward FFT

           ! Determine scaling factor, qfac = exp[iq(r-r')], for current q point
           ! FFTW only does the FFT on G so we must factor out q such that
           ! W_rr'(q) = exp[iq(r-r')]{exp(iGr).W_GG'(q).exp(-iG'r')}
           write(6,'(a8,a15,a15,i4,a3,i4,a4)') "Scaling",adjustr(calcname),"for unfolded q",iqs,"of",nqstar,"..."

           ! Convert q to Cartesian 2pi/a units
           ! Note that qstar is not used for chi since V isn't computed in that case
           qcart = 0.d0
           qcart = qstar(1,iqs)*b(:,1) + qstar(2,iqs)*b(:,2) + qstar(3,iqs)*b(:,3)

           ! Convert q to Cartesian 1/au units
           qcart = qcart * blat

           ! Loop over r and r' in order to multiply by the required exp[q(r-r')]
           ! factor
           ! Nic: supercell modifications
           ! nr -> nr_sc
           ! rcart -> sc_rcart
           ! fftbox -> sc_fftbox
           ! rvec2real -> scvec2real
           call get_sc_fftbox(fftbox,nr_sc,scvec2real,rvec2real,sc2uc_idx,sc_fftbox)
           !$omp parallel do default(private) shared(sc_rcart,scvec2real,qcart,sc_fftbox,nr_sc)
           do ir = 1,nr_sc ! r loop
              do irp = 1,nr_sc ! r' loop
                 qfac = (0.d0,0.d0)
                 r = sc_rcart(:,ir) ! r
                 rp = sc_rcart(:,irp) ! r'
                 rtot = r - rp
                 rq = dot_product(qcart,rtot)
                 qfac = exp(I_D*rq)
                 ! Now scale W(r,r') by the corresponding factor for the current q
                 sc_fftbox(scvec2real(1,ir),scvec2real(2,ir),scvec2real(3,ir),scvec2real(1,irp),scvec2real(2,irp),scvec2real(3,irp)) = &
                      sc_fftbox(scvec2real(1,ir),scvec2real(2,ir),scvec2real(3,ir),scvec2real(1,irp),scvec2real(2,irp),scvec2real(3,irp))  * &
                      qfac
              enddo
           enddo
           !$omp end parallel do
           ! Add W_rr'(q) to the previous value of W_rr'
           ! Nic: fftbox -> sc_fftbox
           data_rrp = data_rrp + sc_fftbox

        enddo ! end qstar loop
     else
        ! We are calculating chi; read chimat file

        ! Read header info
        if (initial_read) then
           open(09,file='chimat',form='unformatted',status='old',iostat=ioerr)
           read(09)
           read(09)
           read(09)
           read(09)
           read(09)
           read(09)
           read(09)
           read(09)
           read(09)
           read(09)
           initial_read = .false.
        endif
        ! Read q-dependent info
        read(09)
        read(09) nmtx
        allocate (matrix_ggp (nmtx,nmtx), stat=allocerr)
        allocate (data_ggp (nmtx,nmtx), stat=allocerr)
        allocate (data_ggp1d (nmtx*nmtx), stat=allocerr)
        do igp = 1,nmtx
           ! read the chimat matrix column-by-column
           read(09) (matrix_ggp(ig,igp),ig=1,nmtx)
        enddo

        ! Calculate specified property
        if (calctype .eq. 4) then
           write(6,*) "Calculating chi(G,G')..."
           ! return chi
           data_ggp = matrix_ggp
        else
           ! There is an error
        endif

        ! Put W_GG(q) data into FFT box in preparation for FFTW
        ! Store W_GG(q) matrix in 1D array suitable for FFTW.
        ! Stored as: [G1G'1, G1G'2, ..., G1G'N, G2G'1, G2G'2, ..., G2G'N, ...,
        ! GNG'1, GNG'2, ..., GNG'N]
        i = 1
        do ig = 1,nmtx
           do igp = 1,nmtx
              data_ggp1d(i) = data_ggp(ig,igp)
              i = i + 1
           enddo
        enddo

        ! Put data into FFT box
        ! Note that this routine places the data assuming G -> G and G' -> -G'
        ! since we have W_rr'(q) = exp[i(q+G)r].W_GG'(q).exp[-i(q+G')r']
        call put_into_fftbox(data_ggp1d, nmtx, ng, glist, fftbox, maxnfft)

        ! Do the backward FFT
        call do_fft(fftbox,maxnfft,1) ! -1 for forward FFT and 1 for backward FFT

        ! Determine scaling factor, qfac = exp[iq(r-r')], for current q point
        ! FFTW only does the FFT on G so we must factor out q such that
        ! W_rr'(q) = exp[iq(r-r')]{exp(iGr).W_GG'(q).exp(-iG'r')}
        write(6,'(a8,a15,a10)') "Scaling",calcname,"for q0..."

        ! Convert q to Cartesian 2pi/a units
        ! Note that qstar is not used for chi since V isn't computed in that case
        qcart = 0.d0
        qcart = q(1,iq)*b(:,1) + q(2,iq)*b(:,2) + q(3,iq)*b(:,3)

        ! Convert q to Cartesian 1/au units
        qcart = qcart * blat

        ! Loop over r and r' in order to multiply by the required exp[q(r-r')]
        ! factor
        !$omp parallel do default(shared) private(qfac,ir,irp,r,rp,rtot,rq)
        do ir = 1,nr_sc ! r loop
           do irp = 1,nr_sc ! r' loop
              qfac = (0.d0,0.d0)
              r = sc_rcart(:,ir) ! r
              rp = sc_rcart(:,irp) ! r'
              rtot = r - rp
              rq = dot_product(qcart,rtot)
              qfac = exp(I_D*rq)

              ! Now scale W(r,r') by the corresponding factor for the current q
              sc_fftbox(scvec2real(1,ir),scvec2real(2,ir),scvec2real(3,ir),scvec2real(1,irp),scvec2real(2,irp),scvec2real(3,irp)) = &
                   sc_fftbox(scvec2real(1,ir),scvec2real(2,ir),scvec2real(3,ir),scvec2real(1,irp),scvec2real(2,irp),scvec2real(3,irp))  * &
                   qfac
           enddo
        enddo
        !$omp end parallel do
        data_rrp = data_rrp + sc_fftbox

        ! Add W_rr'(q) to the previous value of W_rr'
        data_rrp = data_rrp + sc_fftbox
     endif

     write(6,*) repeat("=",79)

     ! Clean up
     if (calctype .ne. 4) then
        deallocate(ekin)
        deallocate(v_gp)
        deallocate(qstar)
        deallocate(isortg)
        deallocate(isortgi)
        !Nic
        deallocate(vinv_gp)
     endif
     deallocate(matrix_ggp)
     !x
     deallocate(epstemp)
     deallocate(data_ggp)
     deallocate(data_ggp1d)
  enddo ! end q loop
  !x
  deallocate(ind)
  deallocate(ph)
  deallocate(rqin)
  deallocate(epsqpt)
  deallocate(indrq)
  deallocate(isortg_fzeps)
  deallocate(q_fzeps)
  !Nic
  deallocate (sc2uc_idx)
  deallocate (scvec2real)
  deallocate (fftbox)
  deallocate (sc_fftbox)
  ! Multiply by 1/Omega = 1/(nq_fullBZ * vcell)
  data_rrp = (1/(nqfbz*vcell)) * data_rrp
  ! Close open files
  close(11)
  close(85)
  close(19)
  close(18)

  !==============================================================================
  !WRITE OUTPUT FILES
  ! KN: Presently sections need to be commented/uncommented based on the desired
  !     output type
  !==============================================================================
  ! BINARY OUTPUT FILE
  ! Nic
  ! Write data_rrp to a binary for analysis
  ! Refer to companion utility read_datarrp.f90
  ! HELP: data_rrp can be very huge for a big supercell
  ! HELP: anything other information to put inside?

  write(6,*) "Writing data_rrp to file"
  open(unit=89,file="data_rrp.bin",form='unformatted',status='replace',iostat=ioerr)

  ! KN: the binary files shouldn't have headers
  ! Write headers
  write(89) "Type of calculation"
  write(89) calcname
  write(89) "Lattice vectors in direct coordinates and alat in bohr"
  do i=1,3
     write(89) (a(i,j),j=1,3)
  enddo
  write(89) alat
  write(89) "Size of 6D FFT grid"
  write(89) (maxnfft(i),i=1,6)
  write(89) "Size of supercell"
  write(89) (sc_size(i),i=1,3)
  write(89)

  ! Save data
  write(89) data_rrp
  close(89)
  !------------------------------
  ! PLANAR AVERAGE OF f(r,r) ALONG z
  ! write(6,*) "Writing output files..."
  ! open(unit=97,file="wrrp.cplx.out",form='formatted',status='replace',iostat=ioerr)
  ! open(unit=98,file="wrrp.real.out",form='formatted',status='replace',iostat=ioerr)
  ! open(unit=99,file="wrrp.abs.out",form='formatted',status='replace',iostat=ioerr)
  ! allocate (planavg (maxnfft(3)), STAT=allocerr)
  ! do z = 1,sc_size(3)*maxnfft(3)
  !   tmpsum = 0.d0
  !   do y = 1,sc_size(2)*maxnfft(2)
  !     do x = 1,sc_size(1)*maxnfft(1)
  !       tmpsum = tmpsum + data_rrp(x,y,z,x,y,z)
  !     enddo
  !   enddo
  !   planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
  !   write(97,*) z, planavg(z)
  !   write(98,*) z, real(planavg(z))
  !   write(99,*) z, abs(planavg(z))
  ! enddo
  ! deallocate (planavg)
  !------------------------------
  ! 2D CONTOUR PLOT OF f(r,r') WITH FIXED r'
  ! KN: Note that r' must be specified in terms of FFT index
  ! write(6,*) "Writing output files..."
  ! open(unit=97,file="wrrp.cplx.out",form='formatted',status='replace',iostat=ioerr)
  ! open(unit=98,file="wrrp.real.out",form='formatted',status='replace',iostat=ioerr)
  ! open(unit=99,file="wrrp.abs.out",form='formatted',status='replace',iostat=ioerr)
  !
  ! do z = 1,maxnfft(3)
  !   do x = 1,maxnfft(1)
  !     !fix electron at bonding site (centre of Si-Si bond)
  !     !get W(r,r') for 110 plane, which cuts xy plane diagonally
  !
  !     ! for 12x12x12 FFT grid
  !     write(97,*) x, z, data_rrp(x,x,z,3,3,3)
  !     write(98,*) x, z, real(data_rrp(x,x,z,3,3,3))
  !     write(99,*) x, z, abs(data_rrp(x,x,z,3,3,3))
  !
  !   enddo
  !   ! write blank lines required for gnuplot contour plot input
  !   write(97,*)
  !   write(98,*)
  !   write(99,*)
  ! enddo

  !==============================================================================
  ! CLEANUP
  deallocate(data_rrp)
  deallocate (sc_rcart)
  write(6,*) "Finished!"

end program wrrp

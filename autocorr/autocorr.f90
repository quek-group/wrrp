!=========================================================================
! Performs the 2D autocorrelation of a wrrp.x data_rrp.bin file
! Initial version: Keian Noori 14/01/2019. 
!=========================================================================

#include "f_defs.h"

program autocorr

  use params_m
  use input_m
  use supercell_m
  use interp_m
  implicit none

!------------------------------------------------------------------------------
  ! Variable declarations
  ! General
  integer :: allocerr,ioerr
  character :: calctype*12,date*8, time*10
  ! For wrrp.x data
  parameter cm2bohr = 188971616.463207
  integer :: ngrid(3),nfft(3),sc_size(3)
  real(DP) :: alat,area,conc
  real(DP) :: delawrrp(3,3),delasc(3,3),awrrp(3,3),auc(3,3),asc(3,3),crossproduct(3)
  complex(DPC), allocatable :: datawrrp(:,:,:)
  ! For nearest neighbours and interpolation
  integer :: inn
  integer :: nn_shift(3,8,1)
  real(DP) :: a(3,3),b_shift(3,1)
  integer, allocatable :: nn_ref(:,:,:)
  real(DP), allocatable :: b_ref(:,:)
  complex(DPC) :: tmpdata,datainterp_shift
  complex(DPC) :: datann_shift(8)
  complex(DPC), allocatable :: datainterp_ref(:),datann_ref(:,:)
  ! For autocorrelation
  integer :: i,j,ir,iphi,igrid
  integer :: ndivar1,ndivar2,ndivr,ndivphi,nextrachg
  integer :: chgpos(3)
  real(DP) :: rmax,r,delr,phi,delphi,zcoor,d,v_ref,v_shift,epsconst,eta
  real(DP) :: centresc(3),width(3),height(3),origin(3),delar1(3),delar2(3),shift(2),shiftpoint(3),chgposbohr(3),z(3),z0(3)
  !real(DP), allocatable :: point(:,:)
  real(DP), allocatable :: refpoint(:,:),vgrid2d(:,:),v(:),vmod(:)
  complex(DPC) :: spacesum,phisum
  complex(DPC), allocatable :: c(:)

!------------------------------------------------------------------------------
  write(6,*) repeat("=",79)
  call date_and_time(DATE=date,TIME=time)
  write(6,*) "autocorr.x -- 2D autocorrelation of wrrp.x output..."
  write(6,*) "version: beta"
  write(6,*) "start time: ",date," at ",time
#ifdef CPLX
  write(6,*) "We are running in COMPLEX mode."
#else 
  write(6,*) "We are running in REAL mode."
#endif 
  write(6,*) repeat("=",79)

!------------------------------------------------------------------------------
! wrrp.x binary data file
  write(6,*) "Reading information from wrrp.x output file..."
  call read_wrrp3d(calctype,alat,awrrp,nfft,sc_size,datawrrp)
  write(6,*) "Data is from a wrrp.x ",calctype,"calculation."
  write(6,'(1x,a,3i4)') "The FFT grid size is ",nfft
  ! Compute unit cell lattice vectors
  auc = alat * awrrp
  write(6,*) "The unit cell lattice vectors (bohr) are:"
  write(6,'(3f12.6)') auc
  ! Compute supercell lattice vectors
  asc(:,1) = sc_size(1) * auc(:,1)
  asc(:,2) = sc_size(2) * auc(:,2)
  asc(:,3) = sc_size(3) * auc(:,3)
  write(6,'(1x,a,3i5)') "The supercell size is ",sc_size
  write(6,*) "The supercell lattice vectors (bohr) are:"
  write(6,'(3f12.6)') asc
  ! Compute area of the supercell
  crossproduct(1) = asc(2,1)*asc(3,2) - (asc(3,1)*asc(2,2))
  crossproduct(2) = asc(3,1)*asc(1,2) - (asc(1,1)*asc(3,2))
  crossproduct(3) = asc(1,1)*asc(2,2) - (asc(2,1)*asc(1,2))
  area = norm2(crossproduct)
  write(6,'(1x,a,f12.4,a)') "The supercell area is:",area," bohr^2."
  ! Deltas of the uc and sc lattice vectors (in bohr)
  delawrrp(:,1) = auc(:,1) / nfft(1)
  delawrrp(:,2) = auc(:,2) / nfft(2)
  delawrrp(:,3) = auc(:,3) / nfft(3)
  delasc(:,1) = asc(:,1) / (sc_size(1)*nfft(1))
  delasc(:,2) = asc(:,2) / (sc_size(2)*nfft(2))
  delasc(:,3) = asc(:,3) / (sc_size(3)*nfft(3))
  write(6,*) repeat("=",79)

!------------------------------------------------------------------------------
! Read "autocorr.inp" input file
  write(6,*) "Reading input file..."
  call read_input_autocorr(zcoor,nextrachg,chgpos)
  write(6,'(1x,a,f10.6,a)') "The potential surface plane is at ",zcoor," bohr."
  write(6,'(1x,a,3i5,a)') "The charge perturbation is located at FFT index",chgpos,"."
  chgposbohr = (chgpos(1)-1)*delasc(:,1) + (chgpos(2)-1)*delasc(:,2) + (chgpos(3)-1)*delasc(:,3)
  write(6,'(1x,a,3f12.6,a)') "The charge perturbation is located at ",chgposbohr," bohr."
  write(6,'(1x,a,i4,a)') "There are ", nextrachg+1, " perturbing (impurity) charges in the supercell."
  conc = ((nextrachg+1)/area) * (cm2bohr**2)
  write(6,'(1x,a,d12.4,a)') "The impurity concentration is :",conc," cm^-2."
  write(6,*) repeat("=",79)

!------------------------------------------------------------------------------
  ! Compute reduced 2D grid in the plane specified by zcoor
  ! KN: We choose the reduced cell (i.e. the portion of the cell over which the
  ! autocorrelation is computed) to be half the area of the full plane. This
  ! allows the maximum r, in C(r), to also be half the width (or height, if
  ! smaller) of the plane without causing r to exceed the cell boundaries.

  write(6,*) "Computing reduced 2D grid for autocorrelation of wrrp.x data..."
  ! Settings for autocorrelation
  ndivar1 = 120 ! # divisions of reduced a1, ar1
  ndivar2 = 120 ! # divisions of reduced a2, ar2
  ndivr = 25
  ndivphi = 20
  centresc = 0.5*asc(:,1) + 0.5*asc(:,2)
  width = 0.5*asc(:,1)
  height = 0.5*asc(:,2)
  origin = centresc - 0.5*width - 0.5*height ! origin of reduced grid
  delar1 = width/ndivar1
  delar2 = height/ndivar2
  rmax = min(norm2(width/2),norm2(height/2))
  delr = rmax/ndivr
  delphi = (2*PI_D)/ndivphi
!  zcoor = 14.172945

  ! Allocations
  allocate(datann_ref(8,ndivar1*ndivar2),stat=allocerr)
  allocate(refpoint(3,ndivar1*ndivar2),stat=allocerr)
  allocate(vgrid2d(3,ndivar1*ndivar2),stat=allocerr)
  allocate(v(ndivar1*ndivar2),stat=allocerr)
  allocate(vmod(ndivar1*ndivar2),stat=allocerr)
  allocate(b_ref(3,ndivar1*ndivar2),stat=allocerr)
  allocate(nn_ref(3,8,ndivar1*ndivar2),stat=allocerr)
  allocate(datainterp_ref(ndivar1*ndivar2),stat=allocerr)
  allocate(c(ndivr),stat=allocerr)

!!TEST: WRITE RAW datawrrp DATA AT SPECIFIED nfft(3) PLANE
  open(91,file="w.raw2d.real.dat",form='formatted',status='replace',iostat=ioerr)
  igrid = 19
  do i = 1,sc_size(1)*nfft(1)
    do j = 1,sc_size(2)*nfft(2)
      write(91,'(i5,i5,f12.6)') i,j,real(datawrrp(i,j,igrid))
      !shiftpoint = (i-1)*delasc(:,1) + (j-1)*delasc(:,2) + (igrid-1)*delasc(:,3)
      !write(6,'(3f12.6)') shiftpoint
    enddo
    write(91,*)
  enddo
  close(91)

  ! Compute the reduced 2D grids
  datann_ref = 0.d0
  igrid = 1
  do i = 1,ndivar1
    do j = 1,ndivar2
      ! Grid for wrrp.x data (at surface of slab)
      refpoint(:,igrid) = origin + (i-1)*delar1 + (j-1)*delar2
      refpoint(3,igrid) = zcoor
      ! Grid for bare Coulomb (at plane of charge perturbation)
      vgrid2d(:,igrid) = refpoint(:,igrid)
      vgrid2d(3,igrid) = chgposbohr(3)
      igrid = igrid + 1
    enddo
  enddo
 
  ! Compute nearest neighbours for all points in the reduced 2D grid
  ! LHS and RHS of Ax = b. These are needed for the determination of the
  ! 8 nearest neighbours. Values are saved to variables a and b since
  ! they're overwritten by the DGESV routine
  a = delawrrp
  b_ref = refpoint
  call nearest_neighbours(a,b_ref,nn_ref)
  ! Interpolate wrrp.x output data to all points in the reduced 2D grid
  igrid = 1
  do i = 1,ndivar1
    do j = 1,ndivar2
      do inn = 1,8
        datann_ref(inn,igrid) = datawrrp(nn_ref(1,inn,igrid),nn_ref(2,inn,igrid),nn_ref(3,inn,igrid))
      enddo
      call tlinterp(delawrrp,refpoint(:,igrid),nn_ref(:,:,igrid),datann_ref(:,igrid),tmpdata)
      datainterp_ref(igrid) = tmpdata
      igrid = igrid + 1
    enddo
  enddo

  ! Write interpolated data to file
  open(99,file="tlinterp2d.real.dat",form='formatted',status='replace',iostat=ioerr)
  igrid = 1
  do i = 1,ndivar1
    do j = 1,ndivar2
      write(99,'(i4,i4,f14.8)') & 
            i,j,real(datainterp_ref(igrid))
      igrid = igrid + 1
    enddo
    write(99,*)
  enddo
  close(99)

  ! Log data
  open(90,file="autocorr.log",form='formatted',status='replace',iostat=ioerr)
  write(90,*) "Unit cell lattice vectors in bohr:"
  write(90,'(3f12.6)'), auc
  write(90,*) "Supercell lattice vectors in bohr:"
  write(90,'(3f12.6)'), asc
  write(90,*) "Centre of supercell:"
  write(90,'(3f12.6)'), centresc
  write(90,*) "Width and height of reduced supercell:"
  write(90,'(3f12.6,3f12.6)'), width, height
  write(90,*) "Origin of reduced supercell:"
  write(90,'(3f12.6)'), origin
  write(90,*) "dela1_reduced"
  write(90,'(3f12.6)'), delar1
  write(90,*) "dela2_reduced"
  write(90,'(3f12.6)'), delar2
  write(90,*) "rmax in bohr"
  write(90,'(f10.6)'), rmax
  write(90,*) "delr in bohr"
  write(90,'(f10.6)'), delr
  write(90,*) "delphi in radians"
  write(90,'(f10.6)'), delphi
  write(90,*) repeat("=",79)
  write(90,'(3f12.6)') (refpoint(:,igrid),igrid = 1,ndivar1*ndivar2)
  write(90,*) repeat("=",79)
  write(90,'(3f12.6)') (vgrid2d(:,igrid),igrid = 1,ndivar1*ndivar2)
  ! Log data

!------------------------------------------------------------------------------
  ! Compute bare Coulomb interaction in the reduced 2D plane of the perturbation
  write(6,*) "Computing bare Coulomb potential of impurities in the reduced 2D grid..."
  !epsconst = 2.45
  epsconst = 1.00
  eta = 0.00
  igrid = 1
  do i = 1,ndivar1
    do j = 1,ndivar2
      !r = norm2(vgrid2d(:,igrid) - chgposbohr)
      r = norm2(refpoint(:,igrid) - chgposbohr)
      v(igrid) = 2/(epsconst*(r+eta))

!write(6,'(a,i3,i3)') "i / j ",i,j
!write(6,'(a,3f12.6)') "current grid point (bohr)",vgrid2d(:,igrid)
!write(6,'(a,3f12.6)') "charge position (bohr)",chgposbohr
!write(6,'(a,f12.6)') "delta_r magnitude (bohr)",r
!write(6,*) "===================="
 
      igrid = igrid + 1
    enddo
  enddo

!  ! Compute model bare Coulomb interaction in the reduced 2D plane of the perturbation
!  write(6,*) "Computing model bare Coulomb potential of impurity distribution in the reduced 2D grid..."
!  igrid = 1
!  do i = 1,ndivar1
!    do j = 1,ndivar2
!      z = vgrid2d(:,igrid)
!      z(3) = vgrid2d(3,igrid)+18.8973 ! d = 1nm ~ 19 bohr
!      z0 = chgposbohr
!      z0(3) = chgposbohr(3) - 5669.18 ! image plane is 300 nm below charge plane
!      r = norm2(z-z0)
!!write(6,'(a,i6)') "igrid: ",igrid
!!write(6,'(a,3f12.6)') "vgrid2d: ",vgrid2d(:,igrid)
!!write(6,'(a,3f12.6)') "z: ",z
!!write(6,'(a,3f12.6)') "z0: ",z0
!!write(6,*) "==============="
!      !vmod(igrid) = ((1-3.9)/(1+3.9))*(27000.2114/(4*r)) ! convert Ha to meV for comparison with Shaffique's paper
!      !vmod(igrid) = 27000.2114/(4*3.9*r) ! convert Ha to meV for comparison with Shaffique's paper
!      d = 5669.18 
!      r = norm2(vgrid2d(:,igrid) - chgposbohr)
!      !vmod(igrid) = 27000.2114/(3.9*(sqrt((r*r) + (18.8973*18.8973)))) ! convert Ha to meV for comparison with Shaffique's paper
!      vmod(igrid) = 27000.2114/(3.9*(sqrt((r*r) + (d*d)))) ! convert Ha to meV for comparison with Shaffique's paper
!
!      igrid = igrid + 1
!    enddo
!  enddo

  ! Write bare Coulomb data to file
  open(97,file="v2d.dat",form='formatted',status='replace',iostat=ioerr)
  !open(94,file="vmod2d.dat",form='formatted',status='replace',iostat=ioerr)
  igrid = 1
  do i = 1,ndivar1
    do j = 1,ndivar2
      write(97,'(i4,i4,f15.6)') i,j,v(igrid)
      !write(94,'(i4,i4,f15.6)') i,j,vmod(igrid)
      igrid = igrid + 1
    enddo
    write(97,*)
    !write(94,*)
  enddo
  close(97)
  !close(94)
  
!------------------------------------------------------------------------------
  ! Compute 2D autocorrelation for wrrp.x output data within reduced grid
  ! C(r) = (1/2pi)*int[0,2pi]{dphi*<<V(r^)*V(0)>>}
  ! See eq 1 in PRB 84 235421 (2011)
  write(6,*) "Computing 2D autocorrelation for wrrp.x data within reduced 2D grid..."
  c = 0.d0
  phisum = 0.d0
  spacesum = 0.d0
  do ir = 1,ndivr
    r = (ir-1)*delr
    do iphi = 1,ndivphi
      phi = (1-iphi)*delphi
      shift(1) = r*cos(phi)
      shift(2) = r*sin(phi)
      igrid = 1
      do i = 1,ndivar1
        do j = 1,ndivar2
          ! Find nearest neighbours of shifted point, V(r^)
          a = delawrrp
          b_shift = 0.d0
          shiftpoint(1) = refpoint(1,igrid) + shift(1)
          shiftpoint(2) = refpoint(2,igrid) + shift(2)
          shiftpoint(3) = refpoint(3,igrid)
          ! Find nearest neigbours of shifted point
          b_shift(:,1) = shiftpoint
          call nearest_neighbours(a,b_shift,nn_shift)  
          ! Get wrrp.x output data at each of the 8 nearest neighbours
          datann_shift = 0.d0
          do inn = 1,8
            datann_shift(inn) = datawrrp(nn_shift(1,inn,1),nn_shift(2,inn,1),nn_shift(3,inn,1))
          enddo

          ! Interpolate data for shiftpoint
          call tlinterp(delawrrp,shiftpoint,nn_shift(:,:,1),datann_shift,tmpdata)
          datainterp_shift = tmpdata
          
          ! Get running total over the reduced spatial 2D grid
          spacesum = spacesum + datainterp_ref(igrid)*datainterp_shift
          
          igrid = igrid +1
        enddo ! j
      enddo  ! i
      ! Average over the reduced spatial 2D grid
      ! KN: SHOULD THIS BE AVERAGED OVER ACTUAL SPACE BY MULTIPLYING BY dA AND
      ! THEN DIVIDING BY TOTAL AREA?

!write(6,*) "spacesumW",abs(spacesum)

      spacesum = spacesum/(ndivar1*ndivar2)

      ! Running total over angle phi
      !phisum = phisum + spacesum
      !phisum = phisum + (delphi*spacesum)
      c(ir) = c(ir) + (delphi*spacesum)

    enddo  ! phi
    ! Radially averaged autocorrelation
    !c(ir) = (1/(2*PI_D)) * phisum
    c(ir) = (1/(2*PI_D)) * c(ir)
  enddo  ! r

  ! Write C(r) to file
  open(98,file="c-r.dat",form='formatted',status='replace',iostat=ioerr)
  do ir = 1,ndivr
    r = (ir-1)*delr
    write(98,'(f12.6,f14.8)'),r,abs(c(ir))
  enddo
  close(98)

  ! Compute 2D autocorrelation, in the reduced grid, for the bare Coulomb potential of the impurity
  ! distribution
  write(6,*) "Computing 2D autocorrelation for bare Coulomb potential within reduced 2D grid..."
  c = 0.d0
  phisum = 0.d0
  spacesum = 0.d0
  do ir = 1,ndivr
    r = (ir-1)*delr
    do iphi = 1,ndivphi
      phi = (1-iphi)*delphi
      shift(1) = r*cos(phi)
      shift(2) = r*sin(phi)
      igrid = 1
      do i = 1,ndivar1
        do j = 1,ndivar2
          !KN: THIS WILL NEED TO HAVE AN EXTRA DIMENSION OF RANK nchg TO ACCOUNT
          !FOR THE ADDITIONAL CHARGES
          !shiftpoint(1) = vgrid2d(1,igrid) + shift(1)
          !shiftpoint(2) = vgrid2d(2,igrid) + shift(2)
          !shiftpoint(3) = vgrid2d(3,igrid)
          shiftpoint(1) = refpoint(1,igrid) + shift(1)
          shiftpoint(2) = refpoint(2,igrid) + shift(2)
          shiftpoint(3) = refpoint(3,igrid)
   
          !v_ref = 2/(norm2(vgrid2d(:,igrid)-chgpos))
          !v_ref = 2/(epsconst*(norm2(refpoint(:,igrid) - chgpos) + eta))
          !v_shift = 2/(epsconst*(norm2(shiftpoint - chgpos) + eta))
          v_ref = 2/(epsconst*(norm2(refpoint(:,igrid) - chgposbohr) + eta))
          v_shift = 2/(epsconst*(norm2(shiftpoint - chgposbohr) + eta))

!write(6,'(a,f12.6)') "shiftV",v_shift
!write(6,'(a,f12.6)') "refV",v_ref

          ! Get running total over the reduced spatial 2D grid
          spacesum = spacesum + v_shift*v_ref
          igrid = igrid +1
        enddo ! j
      enddo  ! i
      ! Average over the reduced spatial 2D grid
      ! KN: SHOULD THIS BE AVERAGED OVER ACTUAL SPACE BY MULTIPLYING BY dA AND
      ! THEN DIVIDING BY TOTAL AREA?

!write(6,*) "spacesumV",abs(spacesum)

      spacesum = spacesum/(ndivar1*ndivar2)

      ! Running total over angle phi
      c(ir) = c(ir) + (delphi*spacesum)

    enddo  ! phi
    ! Radially averaged autocorrelation
    c(ir) = (1/(2*PI_D)) * c(ir)
  enddo  ! r

  ! Write C(r) for bare Coulomb to file
  open(92,file="c-r.bare.dat",form='formatted',status='replace',iostat=ioerr)
  do ir = 1,ndivr
    r = (ir-1)*delr
    write(92,'(f12.6,f14.8)'),r,abs(c(ir))
  enddo
  close(92)

!------------------------------------------------------------------------------
  ! End program
  write(6,*) "Finished!"  

  ! Cleanup
  close(90)
  deallocate(datann_ref)
  deallocate(refpoint)
  deallocate(b_ref)
  deallocate(nn_ref)
  deallocate(datainterp_ref)
  deallocate(c)
  deallocate(vgrid2d)
  deallocate(v)
  deallocate(vmod)

end program autocorr

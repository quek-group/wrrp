!=========================================================================
! Computes triliniear interpolation of a point in 3D space 
! Initial version: Keian Noori 09/10/2018 
!=========================================================================

#include "f_defs.h"

module interp_m

  use params_m
  implicit none

  public  :: nearest_neighbours
  public  :: interp6d
  public  :: tlinterp

contains

!------------------------------
! Find the 8 nearest neighbours (NN) in the FFT box for a set of 3D atomic coordinates
! NN are found by solving Ax=b using Lapack
  subroutine nearest_neighbours(a,b,nn)
    integer :: n,nrhs,lda,ipiv,ldb,info ! for Lapack diagonalization
    integer :: iat,i,j
    real(DP), intent(in) :: a(:,:),b(:,:)
    integer, intent(out) :: nn(:,:,:)
    
    n = 3       ! rank of A matrix = 3 since 3D
    lda = n     ! leading degree of A
    ldb = n     ! leading degree of b
    nrhs = size(b,2) ! number of right hand sides, i.e. number of atoms

!write(6,*) a
!write(6,*) 
!do j = 1,3
!  write(6,*) (a(i,j),i=1,3)
!enddo
!write(6,*) 
!write(6,*) b
!write(6,*) 
!write(6,*) nrhs
!write(6,*) 

    ! Solve Ax=b to get atom position in terms of multiples of A (i.e. unit cell
    ! lattice vectors)
    ! KN: note that A,x, and b are, by definition, column vectors. However, A
    ! and b are input here as row vectors, i.e. each lattice vector / atomic
    ! coordinate occupies a row of its respective matrix. To recover b we must
    ! use the expression A' * x' = b, where b in this expression is made up of
    ! nrhs 3-row column vectors
    call dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
    !call dgesv( n, nrhs, transpose(a), lda, ipiv, transpose(b), ldb, info )

!write(6,*) repeat("=",79)
!write(6,*) b
!write(6,*)
!write(6,*) repeat("=",79)
!!write(6,*) transpose(a)*transpose(b)
!write(6,*) "Size of a is ",size(a,1),size(a,2)
!write(6,*) "Size of b is ",size(b,1),size(b,2)
!!write(6,*) matmul(a,b)
!write(6,*) transpose(a)
!write(6,*)
!write(6,*) transpose(b)

    ! Compute the 8 nearest neigbhours for each atom
    do iat = 1,nrhs
      nn(:,1,iat) = (/floor(b(1,iat)),floor(b(2,iat)),floor(b(3,iat))/)
      nn(:,2,iat) = (/floor(b(1,iat)),ceiling(b(2,iat)),floor(b(3,iat))/)
      nn(:,3,iat) = (/floor(b(1,iat)),floor(b(2,iat)),ceiling(b(3,iat))/)
      nn(:,4,iat) = (/floor(b(1,iat)),ceiling(b(2,iat)),ceiling(b(3,iat))/)
      nn(:,5,iat) = (/ceiling(b(1,iat)),floor(b(2,iat)),floor(b(3,iat))/)
      nn(:,6,iat) = (/ceiling(b(1,iat)),ceiling(b(2,iat)),floor(b(3,iat))/)
      nn(:,7,iat) = (/ceiling(b(1,iat)),floor(b(2,iat)),ceiling(b(3,iat))/)
      nn(:,8,iat) = (/ceiling(b(1,iat)),ceiling(b(2,iat)),ceiling(b(3,iat))/)
  
!write(6,*) repeat("=",79)
!write(6,*) nn(:,:,iat)
    enddo

  end subroutine nearest_neighbours
!------------------------------

! Computes the trilinear intepolation of an intermediate point from the
! surrounding rectangle of known points
  subroutine interp6d(datawrrp,dela,r,rp,nn,nnp,datainterpf)
    integer, intent(in) :: nn(:,:),nnp(:,:)             ! nearest neighbour coordinates of r and r'
    real(DP), intent(in) :: dela(:,:),r(:),rp(:)        ! FFT bin, r and r' coordinates
    complex(DPC), intent(in) :: datawrrp(:,:,:,:,:,:)   ! 6D data from wrrp.x
    complex(DPC), intent(out) :: datainterpf            ! final interpolated data point
    integer :: inn,innp
    complex(DPC) :: datann(8),datainterpi(8),tmpdata    ! data values for each of the NN vertices
    
    ! KN: In order to interpolate the 6D data we first fix a nearest neighbour
    ! for r and then extract datawrrp for every NN of r'. This results in
    ! 8 points: datawrrp(nn1,nn1'),datawrrp(nn1,nn2'),...,datawrrp(nn1,nn8').
    ! Using these 8 points we can perform a standard trilinear interpolation to
    ! get datawrrp(nn1,r'). We repeat this procedure for the remaning NN of r,
    ! obtaining datawrrp(nn2,r'),...,datawrrp(nn8,r'). From *these* 8 points we
    ! can interpolate again to finally get datawrrp(r,r')
    ! KN: Note that we add 1 to each NN index since the origin of datawrrp
    ! starts at index 1
    ! KN: NOTE: THERE MIGHT BE ISSUES AT THE FAR EDGE OF THE CELL; MIGHT NEED TO
    ! USE PERIODIC BOUNDARY CONDITIONS TO COPY ORIGIN VALUES TO FAR EDGE 

!write(6,*) "-----INSIDE-----"
!write(6,'(a,3f6.3)') "r is ",r
!write(6,'(a,3f6.3)') "r' is ",rp

    do inn = 1,8
      datann = 0.d0
      ! Extract datawrrp(nn,nn') for all nn'
      do innp = 1,8
        datann(innp) = datawrrp(nn(1,inn)+1,nn(2,inn)+1,nn(3,inn)+1, &
                           nnp(1,innp)+1,nnp(2,innp)+1,nnp(3,innp)+1)

!write(6,*) "inn/innp:",inn,innp
!write(6,*) datann(innp)

      enddo

      ! Interpolate to get datawrrp(nn1,r')
      call tlinterp(dela,rp,nnp,datann,tmpdata)
      datainterpi(inn) = tmpdata
    enddo

!write(6,*) "-----"
!write(6,*) datainterpi

    ! Interpolate to get datawrrp(r,r')
    call tlinterp(dela,r,nn,datainterpi,datainterpf)

!write(6,*) "-----"
!write(6,*) datainterpf
  
  end subroutine interp6d
!------------------------------

! Standard trilinear interpolation for a point in 3D
  subroutine tlinterp(dela,r,nn,datann,c)
    integer, intent(in) :: nn(:,:)
    real(DP), intent(in) :: dela(:,:),r(:)
    complex(DPC), intent (in) :: datann(:)
    ! Trilinear interpolation variables (c.f. Wikipedia)
    !real(DP) :: a10(3),a11(3),a20(3),a21(3),a30(3),a31(3),delr(3)
    real(DP) :: delr(3),a000(3),a010(3),a001(3),a011(3),a100(3),a110(3),a101(3),a111(3)
!real(DP) :: xd,yd,zd,x0,x1,y0,y1,z0,z1!,dela1,dela2,dela3
    real(DP) :: dela1,dela2,dela3
    complex(DPC) :: c00,c01,c10,c11,c0,c1
    complex(DPC), intent(out) :: c

    ! Find the points above and below x,y,z
    ! Mappings between Wikipedia pts and nn
    ! c000 -> nn(:,1)
    ! c010 -> nn(:,2)
    ! c001 -> nn(:,3)
    ! c011 -> nn(:,4)
    ! c100 -> nn(:,5)
    ! c110 -> nn(:,6)
    ! c101 -> nn(:,7)
    ! c111 -> nn(:,8)

    ! Origin of the NN unit cell in bohr
    a000 = dela(:,1)*nn(1,1) + dela(:,2)*nn(2,1) + dela(:,3)*nn(3,1)
    ! Compute remaining NN points (for debugging purposes)
    a010 = dela(:,1)*nn(1,2) + dela(:,2)*nn(2,2) + dela(:,3)*nn(3,2)
    a001 = dela(:,1)*nn(1,3) + dela(:,2)*nn(2,3) + dela(:,3)*nn(3,3)
    a011 = dela(:,1)*nn(1,4) + dela(:,2)*nn(2,4) + dela(:,3)*nn(3,4)
    a100 = dela(:,1)*nn(1,5) + dela(:,2)*nn(2,5) + dela(:,3)*nn(3,5)
    a110 = dela(:,1)*nn(1,6) + dela(:,2)*nn(2,6) + dela(:,3)*nn(3,6)
    a101 = dela(:,1)*nn(1,7) + dela(:,2)*nn(2,7) + dela(:,3)*nn(3,7)
    a111 = dela(:,1)*nn(1,8) + dela(:,2)*nn(2,8) + dela(:,3)*nn(3,8)
!write(6,'(a,3f6.3)') "a000: ",a000
!write(6,'(a,3f6.3)') "a010: ",a010
!write(6,'(a,3f6.3)') "a001: ",a001
!write(6,'(a,3f6.3)') "a011: ",a011
!write(6,'(a,3f6.3)') "a100: ",a100
!write(6,'(a,3f6.3)') "a110: ",a110
!write(6,'(a,3f6.3)') "a101: ",a101
!write(6,'(a,3f6.3)') "a111: ",a111
    ! Position of r relative to origin
    delr = r - a000
    ! Projection of r along each NN unit cell vector 
    ! KN: Note these correspond to xd, yd, and zd from Wikipedia
    dela1 = dot_product(delr,dela(:,1))
    dela2 = dot_product(delr,dela(:,2))
    dela3 = dot_product(delr,dela(:,3))
!write(6,'(3f6.3)') delr
!write(6,'(f6.3)') dela1
!write(6,'(f6.3)') dela2
!write(6,'(f6.3)') dela3

!write(6,*) "nn: "
!write(6,*) nn
!write(6,*) "a:"
!write(6,*) dela

!    ! Compute the points "below" and "above" r along a1, a2, and a3
!    ! KN: Note that while we name variables x, y, and z, they actually
!    ! correspond to a1, a2, and a3, respectively
!    a10 = dela(:,1)*nn(1,1) + dela(:,2)*nn(2,1) + dela(:,3)*nn(3,1)
!    a11 = dela(:,1)*nn(1,5) + dela(:,2)*nn(2,5) + dela(:,3)*nn(3,5)
!    a20 = a10
!    a21 = dela(:,1)*nn(1,2) + dela(:,2)*nn(2,2) + dela(:,3)*nn(3,2)
!    a30 = a10
!    a31 = dela(:,1)*nn(1,3) + dela(:,2)*nn(2,3) + dela(:,3)*nn(3,3)
!write(6,*) "TEST:"
!write(6,'(3f6.3)') a10
!write(6,'(3f6.3)') a11
!write(6,'(3f6.3)') a20
!write(6,'(3f6.3)') a21
!write(6,'(3f6.3)') a30
!write(6,'(3f6.3)') a31
!write(6,*) "-----"
!    x0 = nn(1,1)*dela(1,1) + nn(2,1)*dela(1,2) + nn(3,1)*dela(1,3)
!    x1 = nn(1,8)*dela(1,1) + nn(2,8)*dela(1,2) + nn(3,8)*dela(1,3)
!    y0 = nn(1,1)*dela(2,1) + nn(2,1)*dela(2,2) + nn(3,1)*dela(2,3)
!    y1 = nn(1,8)*dela(2,1) + nn(2,8)*dela(2,2) + nn(3,8)*dela(2,3)
!    z0 = nn(1,1)*dela(3,1) + nn(2,1)*dela(3,2) + nn(3,1)*dela(3,3)
!    z1 = nn(1,8)*dela(3,1) + nn(2,8)*dela(3,2) + nn(3,8)*dela(3,3)
!
!write(6,*) "NN in bohr:"
!write(6,'(6f12.6)') x0,y0,z0,x1,y1,z1
!write(6,*) "---"
!
!    ! Compute coordinate differences (i.e. relative position between points
!    ! below and above r).
!    ! xd = (x-x0)/(x1-x0), etc
!    xd = (r(1)-x0)/(x1-x0)
!    yd = (r(2)-y0)/(y1-y0)
!    zd = (r(3)-z0)/(z1-z0)
!
!write(6,*) "Derivs:"
!write(6,'(6f8.3)') xd,yd,zd

    ! Interpolate the data
    !c00 = datann(1)*(1-xd) + datann(5)*xd
    !c01 = datann(3)*(1-xd) + datann(7)*xd
    !c10 = datann(2)*(1-xd) + datann(6)*xd
    !c11 = datann(4)*(1-xd) + datann(8)*xd
    c00 = datann(1)*(1-dela1) + datann(5)*dela1
    c01 = datann(3)*(1-dela1) + datann(7)*dela1
    c10 = datann(2)*(1-dela1) + datann(6)*dela1
    c11 = datann(4)*(1-dela1) + datann(8)*dela1
    
    !c0 = c00*(1-yd) + c10*yd 
    !c1 = c01*(1-yd) + c11*yd 
    c0 = c00*(1-dela2) + c10*dela2
    c1 = c01*(1-dela2) + c11*dela2

    !c = c0*(1-zd) + c1*zd
    c = c0*(1-dela3) + c1*dela3

  end subroutine tlinterp

end module interp_m

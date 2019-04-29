!================================================================================================
!   File: supercell.f90 by Cheng Lin Quan Nicholas 
!   Date: 1 March 2018
!================================================================================================

module supercell_m
    
    use params_m
    use mergesort_m
    implicit none
    public  :: get_sc_grid
    public  :: get_sc_fftbox
    private :: boxcheck

contains
!=============================================================================
!   HELP ME: Logic is too complicated?
!=============================================================================
subroutine get_sc_grid(a,nfft,rcart,sc_size,sc_a,sc2uc_idx,scvec2real,sc_rcart)

    integer, intent(in) :: nfft(6), sc_size(3)
    real(DP), intent(in) :: a(3,3) 
    real(DP), intent(in) :: rcart(:,:)
    real(DP), intent(out) :: sc_a(3,3)
    integer, intent(out) :: sc2uc_idx(:), scvec2real(:,:)
    real(DP), intent(out) :: sc_rcart(:,:) !(3,nfft(1)*nfft(2)*nfft(3)*sc_size)

!   To explain mysterious variable names
!   HELP ME: maybe we can name the variables better? 
!   counter variables, error checking
    integer :: ioerr, i, j, k, x, y, z, ii, n, nr, allocerr
    integer :: shift_mul(3), pts, nmax, sc_nfft(3) 
    real(DP) :: point(3),space(3)
    logical :: flag
    real(DP), allocatable :: grid(:,:), sc_grid(:,:)

    write(6,*) "Entering get_sc_grid..."

    ! Number of points in unit cell and sc respectively
    pts = nfft(1)*nfft(2)*nfft(3)
    nmax = sc_size(1)*sc_size(2)*sc_size(3)*pts
    sc_nfft(1) = nfft(1)*sc_size(1)
    sc_nfft(2) = nfft(2)*sc_size(2)
    sc_nfft(3) = nfft(3)*sc_size(3)
    allocate(grid(3,pts), stat=allocerr)
    allocate(sc_grid(3,nmax), stat=allocerr)

    ii = 1
    do i=0,(nfft(1)-1)
        do j=0,(nfft(2)-1)
            do k=0,(nfft(3)-1)
                    grid(1,ii) = dble(i)/nfft(1)
                    grid(2,ii) = dble(j)/nfft(2)
                    grid(3,ii) = dble(k)/nfft(3)
                    ii = ii + 1
            enddo
        enddo
    enddo
     
    !!! Check whether grid multiplied by alat gives correct ans
    ! Compute sc lattice vectors
    do i = 1,3
        sc_a(:,i) = sc_size(i)*a(:,i)
    enddo

    do i=1,pts
        do j=1,3
            grid(j,i) = (grid(j,i) - 0.5)/sc_size(j)
        enddo
    enddo

    ! Shift by unit cell lattice vector
    do i=1,3
        shift_mul(i) = -sc_size(i)/2
        space(i) = 1.0d0/sc_size(i)
    enddo

    n = 1
    do i=1,(-2*shift_mul(1)+1)
        shift_mul(2) = -sc_size(2)/2
        do j=1,(-2*shift_mul(2)+1)
            shift_mul(3) = -sc_size(3)/2
            do k=1,(-2*shift_mul(3)+1)
                do nr=1,pts
                    call boxcheck(shift_mul,nr,space,grid,point,flag)
                    ! If point is in sc, save the index
                    if (flag .eqv. .true.) then
                        sc_grid(1,n) = point(1) 
                        sc_grid(2,n) = point(2)
                        sc_grid(3,n) = point(3)
                        sc2uc_idx(n) = nr
                        
                        ! We need this for the mapping
                        scvec2real(1,n)= mod((n-1)/(sc_nfft(2)*sc_nfft(3))&
                                          ,sc_nfft(1)) + 1
                        scvec2real(2,n)= mod((n-1)/sc_nfft(3),sc_nfft(2)) + 1
                        scvec2real(3,n)= mod(n-1,sc_nfft(3)) + 1 !! Check logic
                        n = n + 1
                    endif
                enddo
                shift_mul(3) = shift_mul(3) + 1
             enddo
             shift_mul(2) = shift_mul(2) + 1
        enddo
        shift_mul(1) = shift_mul(1) + 1
    enddo

    ! To be sure that all the points that are supposed to be in the supercell
    ! are in the supercell...
    if ((n-1) .eq. nmax) then
        print *, "Points ok"
    endif

    print *, "Sorting the data points..."
    ! Sort the data according to ascending x,y,z in direct coordinates
    call mergesort1(sc_grid,sc2uc_idx,nmax,1,nmax,3)    
    call mergesort1(sc_grid,sc2uc_idx,nmax,1,nmax,2)    
    call mergesort1(sc_grid,sc2uc_idx,nmax,1,nmax,1)    
    
    ! For future development...
!    print *, "Interpolating onto coarse grid..."
    ! Interpolation onto a coarse grid in the sc with x,y,z grid points
!    read(11,*,iostat=ioerr) (coarse_ngrid(i),i=1,3)
!    ncoarse=coarse_ngrid(1)*coarse_ngrid(2)*coarse_ngrid(3)
!    allocate(coarse_grid(ncoarse,3))
!    allocate(coarse_data(ncoarse))

!    ii = 1
!    do i=0,(coarse_ngrid(1)-1)
!        do j=0,(coarse_ngrid(2)-1)
!            do k=0,(coarse_ngrid(3)-1)
!                    coarse_grid(ii,1) = real(i)/coarse_ngrid(1)
!                    coarse_grid(ii,2) = real(j)/coarse_ngrid(2)
!                    coarse_grid(ii,3) = real(k)/coarse_ngrid(3)
!                    coarse(:) = coarse_grid(ii,:)
!                    call interpol(coarse,sc_ngrid,nmax,wfnvint,sc_data)
!                    coarse_data(ii) = wfnvint
!                    coarse_grid(ii,1:3) = coarse_grid(ii,1:3) - 0.5
!                    ii = ii + 1
!            enddo
!        enddo
!    enddo
!    
!    ! Write to cube.out file
!    call write_cube(coarse_ngrid,coarse_data,sf_v,sf_origin,na,as,ap,ac)

    ! Compute real space vectors in the supercell
    do i=1,nmax
       sc_rcart(1,i) =sc_grid(1,i)*sc_a(1,1) +sc_grid(2,i)*sc_a(1,2) +&
                      sc_grid(3,i)*sc_a(1,3)                     
       sc_rcart(2,i) =sc_grid(1,i)*sc_a(2,1) +sc_grid(2,i)*sc_a(2,2) +&
                      sc_grid(3,i)*sc_a(2,3)                     
       sc_rcart(3,i) =sc_grid(1,i)*sc_a(3,1) +sc_grid(2,i)*sc_a(3,2) +&
                      sc_grid(3,i)*sc_a(3,3)
    enddo

    ! Clean up
    deallocate(grid)
    deallocate(sc_grid) 
    write(6,*), "Finished get_sc_grid"
end subroutine get_sc_grid
!=============================================================================

!=============================================================================
subroutine get_sc_fftbox(fftbox,nr_sc,scvec2real,rvec2real,sc2uc_idx,sc_fftbox)
    
    complex(DPC), intent(in) :: fftbox(:,:,:,:,:,:)
    integer, intent (in) :: nr_sc, scvec2real(:,:), rvec2real(:,:), &
                            sc2uc_idx(:)
    complex(DPC), intent(in out) :: sc_fftbox(:,:,:,:,:,:)
    integer :: r, rp 

    ! For each data point in sc, map it to the unit cell
!$omp parallel do default(private) shared(sc_fftbox,fftbox,scvec2real,rvec2real,sc2uc_idx,nr_sc)
    do r=1,nr_sc
        do rp=1,nr_sc
            !write(6,*) r, rp, sc2uc_idx(r), sc2uc_idx(rp)
            sc_fftbox(scvec2real(1,r),scvec2real(2,r),scvec2real(3,r),&
                scvec2real(1,rp),scvec2real(2,rp),scvec2real(3,rp)) = &
            fftbox(rvec2real(1,sc2uc_idx(r)),rvec2real(2,sc2uc_idx(r)),rvec2real(3,sc2uc_idx(r)),&
                rvec2real(1,sc2uc_idx(rp)),rvec2real(2,sc2uc_idx(rp)),rvec2real(3,sc2uc_idx(rp)))
        enddo
    enddo
!$omp end parallel do
end subroutine get_sc_fftbox
!=============================================================================

!=============================================================================
!   Check whether a point falls into the supercell
subroutine boxcheck(shift_mul,nr,space,grid,point,flag)
    integer, intent(in) :: shift_mul(3),nr
    real(DP), intent(in) ::  space(3) 
    real(DP), intent(in) :: grid(:,:)
    real(DP), intent(out) :: point(3)
    logical, intent(out) :: flag
    integer :: ii
  
    flag = .true.

    do ii=1,3
        point(ii) = grid(ii,nr) + shift_mul(ii)*space(ii)
    enddo

    ! Supercell defined to be -0.5 < point(i) <= 0.5
    if( point(1) < -0.5 .or. point(1) >= 0.5) then
        flag = .false.
    endif
    if( point(2) < -0.5 .or. point(2) >= 0.5) then
        flag = .false.
    endif
    if( point(3) < -0.5 .or. point(3) >= 0.5) then
        flag = .false.
    endif


end subroutine boxcheck
!=============================================================================

end module supercell_m


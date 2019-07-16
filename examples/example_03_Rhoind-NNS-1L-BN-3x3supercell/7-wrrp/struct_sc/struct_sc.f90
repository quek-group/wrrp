! Parse atomic info and copy to suitable supercell
program struct_sc
    use params_m
    use supercell_m
    implicit none
    integer :: nat, sc_size(3)
    real(DP) :: a(3,3),sc_a(3,3) 
    real(DP) :: point(3)
!------------Look at variables again----------------!
    integer :: ioerr, i, j, k, ii, n, nr, allocerr
    integer :: x, y, z, shift_mul(3), pts, nmax, sc_nfft(3) 
    real(DP) :: space(3)
    logical :: flag
    real(DP), allocatable :: grid(:,:), sc_grid(:,:), sc_rcart(:,:)
    integer, allocatable :: typ(:), sc_typ(:)
    write(6,*) "Entering get_sc_grid..."

    ! Number of points in unit cell
    open(unit=10,file="pts",form="formatted",status="old")
    do i=1,3
        read(10,*) (a(j,i),j=1,3)
    enddo
    read(10,*) nat
    allocate(grid(3,nat))
    allocate(typ(nat))
    do i=1,nat
        read(10,*) typ(i), (grid(j,i),j=1,3)  ! In crystal coords
    enddo
    read(10,*) (sc_size(i),i=1,3)
    
    pts = nat 
    nmax = sc_size(1)*sc_size(2)*sc_size(3)*pts
     
    !!! Check whether grid multiplied by alat gives correct ans
    ! Read supercell (sc) size
    allocate(sc_grid(3,nmax), stat=allocerr)
    allocate(sc_typ (nmax), stat=allocerr)
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
                    if (flag .eqv. .true.) then
                        sc_grid(1,n) = point(1) 
                        sc_grid(2,n) = point(2)
                        sc_grid(3,n) = point(3)
                        sc_typ (n) = typ(nr)
                        n = n + 1
                    endif
                enddo
                shift_mul(3) = shift_mul(3) + 1
             enddo
             shift_mul(2) = shift_mul(2) + 1
        enddo
        shift_mul(1) = shift_mul(1) + 1
    enddo

    if ((n-1) .eq. nmax) then
        print *, "Points ok"
    endif
    call mergesort1(sc_typ,sc_grid,nmax,1,nmax)    
    allocate(sc_rcart(3,nmax))
    do i=1,nmax
       sc_rcart(1,i) =sc_grid(1,i)*sc_a(1,1) +sc_grid(2,i)*sc_a(1,2) +&
                      sc_grid(3,i)*sc_a(1,3)                     
       sc_rcart(2,i) =sc_grid(1,i)*sc_a(2,1) +sc_grid(2,i)*sc_a(2,2) +&
                      sc_grid(3,i)*sc_a(2,3)                     
       sc_rcart(3,i) =sc_grid(1,i)*sc_a(3,1) +sc_grid(2,i)*sc_a(3,2) +&
                      sc_grid(3,i)*sc_a(3,3)
    enddo

    ! Write output
    open(unit=11,file="sc_pts",form="formatted")
    do i=1,3
        write(11,*) (sc_a(j,i),j=1,3)
    enddo
    write(11,*) nmax, " 1"
    do i=1,nmax
        write(11,*) sc_typ(i),(sc_rcart(j,i),j=1,3)
    enddo
    deallocate(grid)
    deallocate(typ)
    deallocate(sc_grid) 
    deallocate(sc_typ)
    deallocate(sc_rcart)
    close(10)
    close(11)
    write(6,*) "Finished get_sc_grid"
end program struct_sc

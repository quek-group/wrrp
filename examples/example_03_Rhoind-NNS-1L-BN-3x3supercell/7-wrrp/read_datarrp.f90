!================================================================================================
!   File: read_datarrp.f90 (v0.1) by Nicholas Lin Quan Cheng
!   Centre for Advanced 2D Materials, National University of Singapore
!   Dept of Physics, Faculty of Science, National University of Singapore 
!   Date: 23 April 2018
!   Read in output binary file from wrrp.x to perform post-processing of data
!================================================================================================
program read_datarrp

    implicit none
    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: DPC = kind((1.0d0,1.0d0))
    character :: calcname*12
    integer :: i, j, k, l, m, n, x, y, z, ioerr, allocerr, maxnfft(6), sc_size(3)
    real(DP) :: alat, a(3,3), a1(3), a2(3), a3(3), zeroes(3), ones
    real(DP), allocatable :: fullchi(:,:,:)
!    complex(DPC) :: tmpsum,tmpsum1,tmpsum2,avg
    complex(DPC), allocatable :: data_rrp(:,:,:,:,:,:), planavg(:)

    ! Read data_rrp.bin
    open(unit=11,file="data_rrp.bin",form='unformatted',status='old',iostat=ioerr)
    
    write(6,*) "This is read_datarrp.cplx.x"
    write(6,*) "Reading data_rrp.bin"
    read(11) ! text
    read(11) calcname
    write(6,*) "calcname:", calcname
    read(11) ! text
    do i=1,3
        read(11) (a(i,j),j=1,3)
    enddo
    read(11) alat
    read(11) ! text
    read(11) (maxnfft(i),i=1,6)
    read(11) ! text
    read(11) (sc_size(i),i=1,3)
    read(11) ! blank line
    allocate (data_rrp (sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3),&
              sc_size(1)*maxnfft(4),sc_size(2)*maxnfft(5),sc_size(3)*maxnfft(6)))

    read(11) data_rrp

    a1 = a(1,:)
    a2 = a(2,:)
    a3 = a(3,:)
    ones = 1.0d0

    write(6,*) "Reading data_rrp.bin complete"
    write(6,*) "For the record, the calculation is for ", calcname
    open(unit=74,file="hbn.belowb.325.xsf",form='formatted',status='replace',iostat=ioerr)
    open(unit=75,file="hbn.belown.325.xsf",form='formatted',status='replace',iostat=ioerr)
    open(unit=76,file="hbn.belowh.325.xsf",form='formatted',status='replace',iostat=ioerr)
    allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)))

    fullchi = 0.d0
    do z = 1,sc_size(3)*maxnfft(3)
      do y = 1,sc_size(2)*maxnfft(2)
        do x = 1,sc_size(1)*maxnfft(1)
          fullchi(x,y,z) = real(data_rrp(x,y,z,16,16,12))
        enddo
      enddo
    enddo
    zeroes = 0.d0
    call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
                         sc_size(3)*maxnfft(3),&
                         ones,ones,ones, &
                         zeroes, a1, a2, a3, alat, 74)

    fullchi = 0.d0
    do z = 1,sc_size(3)*maxnfft(3)
      do y = 1,sc_size(2)*maxnfft(2)
        do x = 1,sc_size(1)*maxnfft(1)
          fullchi(x,y,z) = real(data_rrp(x,y,z,21,21,12))
        enddo
      enddo
    enddo
    zeroes = 0.d0
    call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
                         sc_size(3)*maxnfft(3),&
                         ones,ones,ones, &
                         zeroes, a1, a2, a3, alat, 75)

    fullchi = 0.d0
    do z = 1,sc_size(3)*maxnfft(3)
      do y = 1,sc_size(2)*maxnfft(2)
        do x = 1,sc_size(1)*maxnfft(1)
          fullchi(x,y,z) = real(data_rrp(x,y,z,18,18,12))
        enddo
      enddo
    enddo
    zeroes = 0.d0
    call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
                         sc_size(3)*maxnfft(3),&
                         ones,ones,ones, &
                         zeroes, a1, a2, a3, alat, 76)

    deallocate(fullchi)

    close(74)
    close(75)
    close(76)

    close(11)
    deallocate(data_rrp)
    write(6,*) "Done"
end program read_datarrp

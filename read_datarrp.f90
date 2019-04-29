!================================================================================================
!   File: read_datarrp.f89 (v0.1) by Nicholas Lin Quan Cheng
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
    complex(DPC) :: tmpsum, avg
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

     write(6,*) data_rrp(:,1,1,1,1,1)
  a1 = a(1,:)
  a2 = a(2,:)
  a3 = a(3,:)
  ones = 1.0d0

    write(6,*) "Reading data_rrp.bin complete"
    write(6,*) "For the record, the calculation is for ", calcname
    ! Output files here
  open(unit=97,file="wrrp.cplx.out",form='formatted',status='replace',iostat=ioerr)
  open(unit=98,file="wrrp.real.out",form='formatted',status='replace',iostat=ioerr)
  open(unit=99,file="wrrp.abs.out",form='formatted',status='replace',iostat=ioerr)
  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
  do z = 1,sc_size(3)*maxnfft(3)
    tmpsum = 0.d0
    do y = 1,sc_size(2)*maxnfft(2)
      do x = 1,sc_size(1)*maxnfft(1)
        tmpsum = tmpsum + data_rrp(x,y,z,x,y,z)
      enddo
    enddo
    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
    write(97,*) z, planavg(z)
    write(98,*) z, real(planavg(z))
    write(99,*) z, abs(planavg(z))
  enddo
  deallocate (planavg)

!===========================================================================================
! 6l gr
!  open(unit=109,file="gr.below.plavg.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,5,5,9)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(109,*) z, abs(planavg(z))
!  enddo
!  deallocate (planavg)
!write(6,*) "0a"
!write(6,*) data_rrp(5,5,61,5,5,9)
!write(6,*) "below"
!write(6,*) data_rrp(5,5,17,5,5,9)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,61,x,y,9)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,17,x,y,9)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
!  open(unit=98,file="gr.6l.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  planavg = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,x,y,17)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(98,*) "# top:61, btm:17"
!    write(98,*) z, real(planavg(z))
!  enddo
!===========================================================================================
! 4l gr
!  open(unit=109,file="gr.below.plavg.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,5,5,9)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(109,*) z, abs(planavg(z))
!  enddo
!  deallocate (planavg)
!write(6,*) "0a"
!write(6,*) data_rrp(5,5,41,5,5,9)
!write(6,*) "1a"
!write(6,*) data_rrp(5,5,44,5,5,9)
!write(6,*) "below"
!write(6,*) data_rrp(5,5,16,5,5,9)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,41,x,y,9)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,16,x,y,9)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
!  open(unit=98,file="gr.4l.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  planavg = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,x,y,16)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(98,*) "# top:41, btm:16"
!    write(98,*) z, real(planavg(z))
!  enddo

!===========================================================================================
! 2l gr
!  open(unit=109,file="gr.below.plavg.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,5,5,16)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(109,*) z, abs(planavg(z))
!  enddo
!  deallocate (planavg)
!write(6,*) "0a"
!write(6,*) data_rrp(5,5,32,5,5,16)
!write(6,*) "1a"
!write(6,*) data_rrp(5,5,35,5,5,16)
!write(6,*) "below"
!write(6,*) data_rrp(5,5,23,5,5,16)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,32,x,y,16)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,23,x,y,16)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
!  open(unit=98,file="gr.2l.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  planavg = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,x,y,23)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(98,*) "# top:32, btm:23"
!    write(98,*) z, real(planavg(z))
!  enddo
!===========================================================================================
! 1l gr
!  open(unit=109,file="gr.below.plavg.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,5,5,12)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(109,*) z, abs(planavg(z))
!  enddo
!  deallocate (planavg)
!write(6,*) "0a"
!write(6,*) data_rrp(5,5,20,5,5,12)
!write(6,*) "1a"
!write(6,*) data_rrp(5,5,22,5,5,12)
!write(6,*) "below"
!write(6,*) data_rrp(5,5,18,5,5,12)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,20,x,y,12)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,18,x,y,12)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
  open(unit=98,file="gr.1l.out",form='formatted',status='replace',iostat=ioerr)
  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
  planavg = 0.d0
  do z = 1,sc_size(3)*maxnfft(3)
    tmpsum = 0.d0
    do y = 1,sc_size(2)*maxnfft(2)
      do x = 1,sc_size(1)*maxnfft(1)
        tmpsum = tmpsum + data_rrp(x,y,z,x,y,20)
      enddo
    enddo
    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
    write(98,*) "# top:22, btm:20"
    write(98,*) z, real(planavg(z))
  enddo

!===========================================================================================
! SiO2 41
!  open(unit=71,file="sio2.o.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="sio2.aboveo.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,26,28,39)) ! on C
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,26,28,47)) ! above F
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!
!  open(unit=75,file="sio2.aboveo.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,24,24,47))
!      enddo
!    enddo
!  enddo
!==============================================================================================
! BP
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,38,x,y,24)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,33,x,y,24)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
!  open(unit=71,file="bp.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="bp.below.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,27,15,33)) ! on C
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,27,15,24)) ! on F
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  open(unit=75,file="bp.below.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,24,15,24))
!      enddo
!    enddo
!  enddo
!============================================================================================
! CF 90 lines
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,22,x,y,38)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,30,x,y,38)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
!  open(unit=71,file="cf.c.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="cf.f.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="cf.abovef.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,27)) ! on C
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,30)) ! on F
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,38)) ! above F
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!
!  open(unit=75,file="cf.abovef.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,21,21,38))
!      enddo
!    enddo
!  enddo
!  open(unit=103,file="cf.abovef.24.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="cf.abovef.25.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="cf.abovef.26.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="cf.abovef.27.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="cf.abovef.28.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="cf.abovef.29.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="cf.abovef.30.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="cf.abovef.31.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="cf.abovef.32.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="cf.abovef.33.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,24,17,17,38))
!        write(104,*) x,y, real(data_rrp(x,y,25,17,17,38))
!        write(105,*) x,y, real(data_rrp(x,y,26,17,17,38))
!        write(106,*) x,y, real(data_rrp(x,y,27,17,17,38))
!        write(107,*) x,y, real(data_rrp(x,y,28,17,17,38))
!        write(108,*) x,y, real(data_rrp(x,y,29,17,17,38))
!        write(109,*) x,y, real(data_rrp(x,y,30,17,17,38))
!        write(110,*) x,y, real(data_rrp(x,y,31,17,17,38))
!        write(111,*) x,y, real(data_rrp(x,y,32,17,17,38))
!        write(112,*) x,y, real(data_rrp(x,y,33,17,17,38))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!============================================================================================
! AlN 105 lines
!  open(unit=71,file="aln.al.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="aln.n.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="aln.aboveal.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=74,file="aln.aboven.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,19,21)) ! on Si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,23,25,21)) ! on C
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,19,12)) ! above Ti
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,23,25,12)) ! above O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 74)
!  deallocate(fullchi)
!
!  open(unit=75,file="aln.aboveal.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,21,24,12))
!      enddo
!    enddo
!  enddo
!  open(unit=76,file="aln.al.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(76,*) x, y, z, real(data_rrp(x,y,z,21,24,21))
!      enddo
!    enddo
!  enddo
!  open(unit=103,file="aln.aboveal.16.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="aln.aboveal.17.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="aln.aboveal.18.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="aln.aboveal.19.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="aln.aboveal.20.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="aln.aboveal.21.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="aln.aboveal.22.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="aln.aboveal.23.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="aln.aboveal.24.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="aln.aboveal.25.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,16,17,19,12))
!        write(104,*) x,y, real(data_rrp(x,y,17,17,19,12))
!        write(105,*) x,y, real(data_rrp(x,y,18,17,19,12))
!        write(106,*) x,y, real(data_rrp(x,y,19,17,19,12))
!        write(107,*) x,y, real(data_rrp(x,y,20,17,19,12))
!        write(108,*) x,y, real(data_rrp(x,y,21,17,19,12))
!        write(109,*) x,y, real(data_rrp(x,y,22,17,19,12))
!        write(110,*) x,y, real(data_rrp(x,y,23,17,19,12))
!        write(111,*) x,y, real(data_rrp(x,y,24,17,19,12))
!        write(112,*) x,y, real(data_rrp(x,y,25,17,19,12))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!============================================================================================
! GaN 105 lines

!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,25,x,y,12)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "top"
!    write(6,*) real(avg)
!  avg = 0.d0
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,17,x,y,12)
!      enddo
!    enddo
!    avg = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(6,*) "bottom"
!    write(6,*) real(avg)
!  open(unit=71,file="gan.ga.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="gan.n.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="gan.abovega.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=74,file="gan.aboven.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,19,21)) ! on Si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,23,25,21)) ! on C
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,19,12)) ! above Ti
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,23,25,12)) ! above O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 74)
!  deallocate(fullchi)
!
!  open(unit=75,file="gan.abovega.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,21,24,12))
!      enddo
!    enddo
!  enddo
!  open(unit=103,file="gan.abovega.16.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="gan.abovega.17.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="gan.abovega.18.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="gan.abovega.19.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="gan.abovega.20.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="gan.abovega.21.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="gan.abovega.22.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="gan.abovega.23.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="gan.abovega.24.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="gan.abovega.25.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,16,17,19,12))
!        write(104,*) x,y, real(data_rrp(x,y,17,17,19,12))
!        write(105,*) x,y, real(data_rrp(x,y,18,17,19,12))
!        write(106,*) x,y, real(data_rrp(x,y,19,17,19,12))
!        write(107,*) x,y, real(data_rrp(x,y,20,17,19,12))
!        write(108,*) x,y, real(data_rrp(x,y,21,17,19,12))
!        write(109,*) x,y, real(data_rrp(x,y,22,17,19,12))
!        write(110,*) x,y, real(data_rrp(x,y,23,17,19,12))
!        write(111,*) x,y, real(data_rrp(x,y,24,17,19,12))
!        write(112,*) x,y, real(data_rrp(x,y,25,17,19,12))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!============================================================================================
! SiC 105 lines
!  open(unit=71,file="sic.si.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="sic.c.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="sic.abovesi.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=74,file="sic.abovec.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,19,21)) ! on Si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,23,25,21)) ! on C
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,19,12)) ! above Ti
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,23,25,12)) ! above O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 74)
!  deallocate(fullchi)
!
!  open(unit=75,file="sic.abovesi.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,21,24,12))
!      enddo
!    enddo
!  enddo
!  open(unit=103,file="sic.abovesi.16.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="sic.abovesi.17.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="sic.abovesi.18.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="sic.abovesi.19.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="sic.abovesi.20.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="sic.abovesi.21.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="sic.abovesi.22.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="sic.abovesi.23.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="sic.abovesi.24.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="sic.abovesi.25.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,16,17,19,12))
!        write(104,*) x,y, real(data_rrp(x,y,17,17,19,12))
!        write(105,*) x,y, real(data_rrp(x,y,18,17,19,12))
!        write(106,*) x,y, real(data_rrp(x,y,19,17,19,12))
!        write(107,*) x,y, real(data_rrp(x,y,20,17,19,12))
!        write(108,*) x,y, real(data_rrp(x,y,21,17,19,12))
!        write(109,*) x,y, real(data_rrp(x,y,22,17,19,12))
!        write(110,*) x,y, real(data_rrp(x,y,23,17,19,12))
!        write(111,*) x,y, real(data_rrp(x,y,24,17,19,12))
!        write(112,*) x,y, real(data_rrp(x,y,25,17,19,12))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!============================================================================================
! TiO2 105 lines
!  open(unit=71,file="tio2.ti.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="tio2.o.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="tio2.aboveti.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=74,file="tio2.aboveo.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,9,19,24)) ! on Ti
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,9,32,27)) ! on O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,9,19,32)) ! above Ti
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,9,32,35)) ! above O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 74)
!  deallocate(fullchi)
!
!  open(unit=75,file="tio2.aboveti.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,13,28,35))
!      enddo
!    enddo
!  enddo
!  open(unit=103,file="tio2.aboveti.21.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="tio2.aboveti.22.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="tio2.aboveti.23.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="tio2.aboveti.24.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="tio2.aboveti.25.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="tio2.aboveti.26.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="tio2.aboveti.27.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="tio2.aboveti.28.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="tio2.aboveti.29.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="tio2.aboveti.30.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,21,9,19,32))
!        write(104,*) x,y, real(data_rrp(x,y,22,9,19,32))
!        write(105,*) x,y, real(data_rrp(x,y,23,9,19,32))
!        write(106,*) x,y, real(data_rrp(x,y,24,9,19,32))
!        write(107,*) x,y, real(data_rrp(x,y,25,9,19,32))
!        write(108,*) x,y, real(data_rrp(x,y,26,9,19,32))
!        write(109,*) x,y, real(data_rrp(x,y,27,9,19,32))
!        write(110,*) x,y, real(data_rrp(x,y,28,9,19,32))
!        write(111,*) x,y, real(data_rrp(x,y,29,9,19,32))
!        write(112,*) x,y, real(data_rrp(x,y,30,9,19,32))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!============================================================================================
! Mgo 105 lines
!  open(unit=71,file="mgo.mg.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="mgo.o.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="mgo.abovemg.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=74,file="mgo.aboveo.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,19,19,40)) ! on Mg
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,19,40)) ! on O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,19,19,48)) ! above Mg
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,19,48)) ! above O
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 74)
!  deallocate(fullchi)
!
!  open(unit=75,file="mgo.abovemg.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(75,*) x, y, z, real(data_rrp(x,y,z,19,19,48))
!      enddo
!    enddo
!  enddo
!  open(unit=103,file="mgo.abovemg.35.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="mgo.abovemg.36.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="mgo.abovemg.37.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="mgo.abovemg.38.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="mgo.abovemg.39.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="mgo.abovemg.40.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="mgo.abovemg.41.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="mgo.abovemg.42.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="mgo.abovemg.43.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="mgo.abovemg.44.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,35,19,19,48))
!        write(104,*) x,y, real(data_rrp(x,y,36,19,19,48))
!        write(105,*) x,y, real(data_rrp(x,y,37,19,19,48))
!        write(106,*) x,y, real(data_rrp(x,y,38,19,19,48))
!        write(107,*) x,y, real(data_rrp(x,y,39,19,19,48))
!        write(108,*) x,y, real(data_rrp(x,y,40,19,19,48))
!        write(109,*) x,y, real(data_rrp(x,y,41,19,19,48))
!        write(110,*) x,y, real(data_rrp(x,y,42,19,19,48))
!        write(111,*) x,y, real(data_rrp(x,y,43,19,19,48))
!        write(112,*) x,y, real(data_rrp(x,y,44,19,19,48))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo

!==============================================================================================
! si100 551 96 lines
!  open(unit=71,file="si.si.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="si.abovesi.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,25,25,48)) ! on si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,25,25,56)) ! above si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!    open(unit=74,file="si.si.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!        do y = 1,sc_size(2)*maxnfft(2)
!      write(74,*)  x, y, real(data_rrp(x,y,48,25,25,48)) 
!        enddo
!        write(74,*) 
!    enddo
!
!  open(unit=73,file="rhoind.abovesi.plavg.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,25,25,56)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(73,*) z, real(planavg(z))
!  enddo
!    close(73)
!
!  open(unit=103,file="si.abovesi.44.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="si.abovesi.45.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="si.abovesi.46.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="si.abovesi.47.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="si.abovesi.48.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="si.abovesi.49.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="si.abovesi.50.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="si.abovesi.51.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="si.abovesi.52.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="si.abovesi.53.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,44,25,25,56))
!        write(104,*) x,y, real(data_rrp(x,y,45,25,25,56))
!        write(105,*) x,y, real(data_rrp(x,y,46,25,25,56))
!        write(106,*) x,y, real(data_rrp(x,y,47,25,25,56))
!        write(107,*) x,y, real(data_rrp(x,y,48,25,25,56))
!        write(108,*) x,y, real(data_rrp(x,y,49,25,25,56))
!        write(109,*) x,y, real(data_rrp(x,y,50,25,25,56))
!        write(110,*) x,y, real(data_rrp(x,y,51,25,25,56))
!        write(111,*) x,y, real(data_rrp(x,y,52,25,25,56))
!        write(112,*) x,y, real(data_rrp(x,y,53,25,25,56))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!  open(unit=113,file="si.abovesi.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(113,*) x, y, z, real(data_rrp(x,y,z,25,25,56))
!      enddo
!    enddo
!  enddo
!==============================================================================================
! si100 96 lines
!  open(unit=71,file="si.si.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="si.abovesi.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,13,48)) ! on si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,13,56)) ! above si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!    open(unit=74,file="si.si.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!        do y = 1,sc_size(2)*maxnfft(2)
!      write(74,*)  x, y, real(data_rrp(x,y,48,13,13,48)) 
!        enddo
!        write(74,*) 
!    enddo
!
!  open(unit=73,file="rhoind.abovesi.plavg.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,13,13,56)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(73,*) z, real(planavg(z))
!  enddo
!    close(73)
!
!  open(unit=103,file="si.abovesi.44.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="si.abovesi.45.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="si.abovesi.46.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="si.abovesi.47.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="si.abovesi.48.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="si.abovesi.49.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="si.abovesi.50.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="si.abovesi.51.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="si.abovesi.52.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="si.abovesi.53.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,44,13,13,56))
!        write(104,*) x,y, real(data_rrp(x,y,45,13,13,56))
!        write(105,*) x,y, real(data_rrp(x,y,46,13,13,56))
!        write(106,*) x,y, real(data_rrp(x,y,47,13,13,56))
!        write(107,*) x,y, real(data_rrp(x,y,48,13,13,56))
!        write(108,*) x,y, real(data_rrp(x,y,49,13,13,56))
!        write(109,*) x,y, real(data_rrp(x,y,50,13,13,56))
!        write(110,*) x,y, real(data_rrp(x,y,51,13,13,56))
!        write(111,*) x,y, real(data_rrp(x,y,52,13,13,56))
!        write(112,*) x,y, real(data_rrp(x,y,53,13,13,56))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!  open(unit=113,file="si.abovesi.3d.out",form='formatted',status='replace',iostat=ioerr)
!  do x = 1,sc_size(1)*maxnfft(1)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do z = 1,sc_size(3)*maxnfft(3)
!        write(113,*) x, y, z, real(data_rrp(x,y,z,19,19,56))
!      enddo
!    enddo
!  enddo
!============================================================================================
!  open(unit=71,file="si.si.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="si.abovesi.xsf",form='formatted',status='replace',iostat=ioerr)
!        allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)&
!            ))
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,15,31,46)) ! on si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,15,31,54)) ! above Si
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!    open(unit=72,file="rhoind.si.plavg.out",form='formatted',status='replace',iostat=ioerr)
!    open(unit=73,file="rhoind.abovesi.plavg.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!        do y = 1,sc_size(2)*maxnfft(2)
!      write(72,*)  x, y, real(data_rrp(x,y,46,15,31,46)) 
!      write(73,*)  x, y, real(data_rrp(x,y,46,15,31,54)) 
!        enddo
!        write(72,*) 
!        write(73,*) 
!    enddo
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,15,31,54)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(73,*) z, real(planavg(z))
!  enddo
!    close(73)
!    close(72)
!
!  open(unit=103,file="si.abovesi.43.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="si.abovesi.44.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="si.abovesi.45.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="si.abovesi.46.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="si.abovesi.47.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="si.abovesi.48.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=109,file="si.abovesi.49.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=110,file="si.abovesi.50.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=111,file="si.abovesi.51.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=112,file="si.abovesi.52.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(1)*maxnfft(1)
!      do y = 1,sc_size(2)*maxnfft(2)
!        write(103,*) x,y, real(data_rrp(x,y,43,15,31,54))
!        write(104,*) x,y, real(data_rrp(x,y,44,15,31,54))
!        write(105,*) x,y, real(data_rrp(x,y,45,15,31,54))
!        write(106,*) x,y, real(data_rrp(x,y,46,15,31,54))
!        write(107,*) x,y, real(data_rrp(x,y,47,15,31,54))
!        write(108,*) x,y, real(data_rrp(x,y,48,15,31,54))
!        write(109,*) x,y, real(data_rrp(x,y,49,15,31,54))
!        write(110,*) x,y, real(data_rrp(x,y,50,15,31,54))
!        write(111,*) x,y, real(data_rrp(x,y,51,15,31,54))
!        write(112,*) x,y, real(data_rrp(x,y,52,15,31,54))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!        write(109,*)
!        write(110,*)
!        write(111,*)
!        write(112,*)
!    enddo
!=======================================================================================
! grbn
!  open(unit=71,file="grbn.gr.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=72,file="grbn.abovegr.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=73,file="grbn.grH.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=74,file="grbn.abovegrH.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=75,file="grbn.B.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=76,file="grbn.belowB.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=77,file="grbn.N.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=78,file="grbn.belowN.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=79,file="grbn.H.xsf",form='formatted',status='replace',iostat=ioerr)
!  open(unit=80,file="grbn.belowH.xsf",form='formatted',status='replace',iostat=ioerr)
!
!  allocate (fullchi(sc_size(1)*maxnfft(1),sc_size(2)*maxnfft(2),sc_size(3)*maxnfft(3)))
!    
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,21,35)) ! on Gr
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 71)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,21,43)) ! above Gr
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 72)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,35)) ! on GrH
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 73)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,43)) ! above GrH
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 74)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,21,1)) ! on B
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 75)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,13,21,65)) ! below B
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 76)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,18,16,1)) ! on N
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 77)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,18,16,65)) ! below N 
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 78)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,1)) ! on H
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 79)
!
!  fullchi = 0.d0
!  do z = 1,sc_size(3)*maxnfft(3)
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        fullchi(x,y,z) = real(data_rrp(x,y,z,17,17,65)) ! below H
!      enddo
!    enddo
!  enddo
!  zeroes = 0.d0
!  call xsf_datagrid_3d(fullchi, sc_size(1)*maxnfft(1), sc_size(2)*maxnfft(2),&
!sc_size(3)*maxnfft(3),&
!            ones,ones,ones, &
!            zeroes, a1, a2, a3, alat, 80)
!
!  deallocate(fullchi)
!    open(unit=82,file="rhoind.aboveGr.out",form='formatted',status='replace',iostat=ioerr)
!    open(unit=83,file="rhoind.belowBN.out",form='formatted',status='replace',iostat=ioerr)
!  allocate (planavg (sc_size(3)*maxnfft(3)), STAT=allocerr)
!
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,13,21,43)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(83,*) z, real(planavg(z))
!  enddo
!  do z = 1,sc_size(3)*maxnfft(3)
!    tmpsum = 0.d0
!    do y = 1,sc_size(2)*maxnfft(2)
!      do x = 1,sc_size(1)*maxnfft(1)
!        tmpsum = tmpsum + data_rrp(x,y,z,18,16,65)
!      enddo
!    enddo
!    planavg(z) = tmpsum / (sc_size(1)*sc_size(2)*maxnfft(1)*maxnfft(2))
!    write(82,*) z, real(planavg(z))
!  enddo
!  deallocate (planavg)
!
!  open(unit=103,file="gr.aboveGr.37.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=104,file="gr.aboveGr.36.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=105,file="gr.aboveGr.35.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=106,file="gr.aboveGr.34.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=107,file="gr.aboveGr.33.out",form='formatted',status='replace',iostat=ioerr)
!  open(unit=108,file="gr.aboveGr.32.out",form='formatted',status='replace',iostat=ioerr)
!    do x = 1,sc_size(2)*maxnfft(2)
!      do y = 1,sc_size(1)*maxnfft(1)
!        write(103,*) x,y, real(data_rrp(x,y,37,13,21,43))
!        write(104,*) x,y, real(data_rrp(x,y,36,13,21,43))
!        write(105,*) x,y, real(data_rrp(x,y,35,13,21,43))
!        write(106,*) x,y, real(data_rrp(x,y,34,13,21,43))
!        write(107,*) x,y, real(data_rrp(x,y,33,13,21,43))
!        write(108,*) x,y, real(data_rrp(x,y,32,13,21,43))
!      enddo
!      write(103,*)
!      write(104,*)
!      write(105,*)
!      write(106,*)
!      write(107,*)
!      write(108,*)
!    enddo

    close(11)
    deallocate(data_rrp)
    write(6,*) "Done"
end program read_datarrp


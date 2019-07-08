!=========================================================================
!   File: mergesort.f90
!   Date: March 2018
!   Merge sort routines to be used in supercell.f90
!
! This code was written at the National University of Singapore (NUS)
! v0.2 N Cheng
!=========================================================================

module mergesort_m
  implicit none
  public :: mergesort1
  private :: merge1

contains
  subroutine merge1(lst, lst2, n, a, middle, b,col)
    integer, intent(in) :: n,a,middle,b,col
    integer :: ai,bi,ti,x
    real*8, intent(in out) :: lst(3,n)
    real*8 :: tmp(3,n)
    integer, intent(in out) :: lst2(n)
    integer :: tmp2(n)

    ai = a
    bi = middle
    ti = a

    do while ((ai .lt. middle) .or. (bi .lt. b))
       if (ai .eq. middle) then
          tmp(:,ti+1) = lst(:,bi+1)
          tmp2(ti+1) = lst2(bi+1)
          bi = bi + 1
       else if (bi .eq. b) then
          tmp(:,ti+1) = lst(:,ai+1)
          tmp2(ti+1) = lst2(ai+1)
          ai = ai + 1
       else if (lst(col,ai+1) .le. lst(col,bi+1)) then
          tmp(:,ti+1) = lst(:,ai+1)
          tmp2(ti+1) = lst2(ai+1)
          ai = ai + 1
       else
          tmp(:,ti+1) = lst(:,bi+1)
          tmp2(ti+1) = lst2(bi+1)
          bi = bi + 1
       end if
       ti = ti + 1
    end do
    do x = a, b - 1
       lst(:,x+1) = tmp(:,x+1)
       lst2(x+1)  = tmp2(x+1)
    end do
  end subroutine merge1

  recursive subroutine mergesort1(lst,lst2,n,a,b,col)
    integer, intent(in):: n,a,b,col
    real*8, intent(in out) :: lst(n,3)
    integer, intent(in out) :: lst2(n)
    integer :: diff
    diff = b - a
    if (diff .lt. 2) then
       return
    else
       diff = diff / 2
       call mergesort1(lst, lst2, n, a, a + diff,col)
       call mergesort1(lst, lst2, n, a + diff, b,col)
       call merge1(lst, lst2, n, a, a + diff, b,col)
    endif
  end subroutine mergesort1

end module mergesort_m

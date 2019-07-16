module mergesort_m
    implicit none
    public :: mergesort1
    private :: merge1

contains
subroutine merge1(lst, lst2, n, a, middle, b)
    integer, intent(in) :: n,a,middle,b
    integer :: ai,bi,ti,x
    real*8, intent(in out) :: lst2(3,n)
    real*8 :: tmp(3,n)
    integer, intent(in out) :: lst(n)
    integer :: tmp2(n)
    
    ai = a
    bi = middle
    ti = a

    do while ((ai .lt. middle) .or. (bi .lt. b))
        if (ai .eq. middle) then
            tmp(:,ti+1) = lst2(:,bi+1)
            tmp2(ti+1) = lst(bi+1)
            bi = bi + 1
        else if (bi .eq. b) then
            tmp(:,ti+1) = lst2(:,ai+1)
            tmp2(ti+1) = lst(ai+1)
            ai = ai + 1
        else if (lst(ai+1) .le. lst(bi+1)) then
            tmp(:,ti+1) = lst2(:,ai+1)
            tmp2(ti+1) = lst(ai+1)
            ai = ai + 1
        else
            tmp(:,ti+1) = lst2(:,bi+1)
            tmp2(ti+1) = lst(bi+1)
            bi = bi + 1
        end if
        ti = ti + 1
    end do
    do x = a, b - 1
            lst2(:,x+1) = tmp(:,x+1)
        lst(x+1)  = tmp2(x+1)
    end do
end subroutine merge1

recursive subroutine mergesort1(lst,lst2,n,a,b)
    integer, intent(in):: n,a,b
    real*8, intent(in out) :: lst2(n,3)
    integer, intent(in out) :: lst(n)
    integer :: diff
    diff = b - a
    if (diff .lt. 2) then
    return
    else
    diff = diff / 2
    call mergesort1(lst, lst2, n, a, a + diff)
    call mergesort1(lst, lst2, n, a + diff, b)
    call merge1(lst, lst2, n, a, a + diff, b )
    endif
end subroutine mergesort1

end module mergesort_m

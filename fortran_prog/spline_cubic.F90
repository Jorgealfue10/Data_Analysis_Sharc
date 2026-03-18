module spline_module
    implicit none
    private 
    public :: spline1d, init_spline, eval_spline
    
    type :: spline1d
        real(kind=8), allocatable :: x(:)
        real(kind=8), allocatable :: y(:)
        real(kind=8), allocatable :: y2(:)
    end type spline1d
    
    contains

    !---------------------------------------------------------------------------!
    ! init_spline
    !---------------------------------------------------------------------------!
    subroutine init_spline(s, x_in, y_in)
        implicit none
        type(spline1d), intent(inout) :: s
        real(kind=8), intent(in) :: x_in(:), y_in(:)

        integer :: n

        n = size(x_in)
        allocate(s%x(n), s%y(n), s%y2(n))

        s%x = x_in
        s%y = y_in
        
        call spline(s%x,s%y,n,1.0d30,1.0d30,s%y2)

    end subroutine init_spline

    !---------------------------------------------------------------------------!
    ! eval_spline
    !---------------------------------------------------------------------------!
    subroutine eval_spline(s, xval, yval)
        implicit none
        type(spline1d), intent(in) :: s
        real(kind=8), intent(in) :: xval
        real(kind=8), intent(out) :: yval
        
        call splint(s%x,s%y,s%y2,size(s%x),xval,yval)
    
    end subroutine eval_spline

    !---------------------------------------------------------------------------!
    ! Cubic spline (Numerical Recipes)
    !---------------------------------------------------------------------------!
    subroutine spline(x,y,n,yp1,ypn,y2)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: yp1, ypn, x(n), y(n)
        real(kind=8), intent(out) :: y2(n)

        real(kind=8) :: u(n)
        integer :: i,k
        real(kind=8) :: p,qn,sig,un
        
        if (yp1.gt.0.99d30) then
            y2(1)=0.d0
            u(1)=0.d0
        else
            y2(1)=-0.5d0
            u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        endif

        do i=2,n-1
            sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
            p = sig*y2(i-1)+2.0d0
            y2(i) = (sig-1.0d0)/p
            u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
                /(x(i+1)-x(i-1)) - sig*u(i-1)) / p
        enddo

        if (ypn.gt.0.99d30) then
            qn=0.d0
            un=0.d0
        else
            qn=0.5d0
            un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif

        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
        do k=n-1,1,-1
            y2(k)=y2(k)*y2(k+1)+u(k)
        enddo
    end subroutine spline
    
    !---------------------------------------------------------------------------!
    ! Evaluación spline
    !---------------------------------------------------------------------------!
    subroutine splint(xa,ya,y2a,n,x,y)
        implicit none

        integer, intent(in) :: n
        real(kind=8), intent(in) :: xa(n),ya(n),y2a(n),x
        real(kind=8), intent(out) :: y

        integer :: klo,khi,k
        real(kind=8) :: h,a,b

        if (x <= xa(1)) then
            h = xa(2) - xa(1)
            y = ((ya(2)-ya(1))/h - y2a(1)*h/6.0d0)* (x-xa(1)) + ya(1)
            return
        endif

        if (x >= xa(n)) then
            h = xa(n) - xa(n-1)
            y = ((ya(n)-ya(n-1))/h - y2a(n)*h/6.0d0)* (x-xa(n)) + ya(n)
            return
        endif

        klo = 1
        khi = n

        do while (khi - klo > 1)
            k = (khi + klo)/2
            if (xa(k) > x) then
                khi = k
            else
                klo = k
            endif
        enddo
        
        h = xa(khi) - xa(klo)

        if (h==0.0d0) stop "Bad xa input"
        
        a = (xa(khi)-x)/h
        b = (x-xa(klo))/h
        
        y = a*ya(klo) + b*ya(khi) + &
            ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * h*h/6.0d0
    end subroutine splint
end module spline_module
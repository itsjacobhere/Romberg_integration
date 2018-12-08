!This module establishes parameters, reals, and integers
!solution to this problem is 128.5
module vars
    integer, parameter :: wp = selected_real_kind(15)
    real(wp), parameter :: alpha=0.3204, beta=2.75, lambda=2.96,ll=exp(lambda),epsi =  0.0001!0.000001 
    real(wp) :: valu,ul,percent
    integer :: maxit, n
    !ll=lower limit
    !ul=upper limit, solving for this
    !epsi= exit criteria
    !n=iteration
end module vars
!________________________________________________________________________
module sub 
use vars
implicit none
    contains
    subroutine romberg(ll,ul,f,valu,epsi,maxit,n)
        implicit none
        !Global variables
        !-----------------------------
        real(wp), intent(in) :: ll, epsi
        real(wp), intent(out) :: valu
        real(wp), intent(inout) :: ul
        integer, intent(in) :: maxit
        integer, intent(out) :: n
        
        !interface external function here
        Interface
            function f(xa)
            use vars
            implicit none
            real(wp), intent(in) :: xa
            real(wp) :: f
            end function f
        end interface
        
        !local var
        !-----------------------------
        real(wp), dimension(maxit+1,maxit+1) :: t
        real(wp) :: ratio, step, s1, s2, s2sum, s2i, sfinal 
        integer :: k,i, j, nt
        
        !initialize local vars
        nt = 2
        ratio = 2 !long/short
        t=0
        
        step = (ul-ll)/nt
        s1 = (f(ll) + f(ul))/2
        
        s2 = 0 !intialize s2 before sum
        do k = 1, (nt-1) !Sum s2
            s2sum = f(ll+ (k*step))
            s2 = s2 + s2sum
        end do

        t(1,1) = step*(s1+s2) 
        n=1 !iterations
        do 
            nt = 2*nt !double traps every time
            step = (ul-ll)/nt
            
            sfinal = 0 !initalize
            !summing s2 here before extrapolations
            do i = 1, (nt/2) !this sum is getting stuck
                s2i= f(ll + step*((2.0*i) - 1.0))
                sfinal = sfinal + s2i
                
            end do
        
            s2 = s2 + sfinal
            t(n+1,1) = step*(s1 + s2)
        
            do j = 2, (n+1) !extrapolations formula
                t(n+1,j) = (ratio**(2.0*(j-1))*t(n+1,j-1)-t(n,j-1))/(ratio**(2.0*(j-1))-1)
            end do

            !check each iteration against the stopping criteria
            if (abs(t(n+1,n) - t(n+1, n+1)) < 0.000001) then !Converged
                valu = t(n+1,n+1)
                exit
            else if (n > maxit) then ! too many iterations
                write(*,*) "The max iterations allowed were hit! n = ",n
                exit
            end if
            n=n+1
        end do
    end subroutine romberg
end module sub
!________________________________________________________________________
!This module contains all the input and output

!________________________________________________________________________
! Written by Jacob Turner
! This program will determine the design flow for a disinfection system
! by integrating a PDF over a range of flows
!________________________________________________________________________
program romberg_integration
    use vars
    use sub

    implicit none
    real(wp) :: answer
    real(wp) :: max_flow, perc

    !interface for function here
    !------------------------------
    Interface
        function f(xa)
            use vars
            implicit none
            real(wp), intent(in) :: xa
            real(wp) :: f
        end function f
    end interface

    !ask user for probability percent and an initial guess
    write(*,*)
    write(*,*) "Enter the desired probability percentile"
    read(*,*) answer
    percent = answer/100
    write(*,*)
    write(*,'(A,F4.2)') "Your desired percentile is ", percent
    write(*,*)
    maxit = 30
    write(*,'(A,I2,A)') "Stopping after ", maxit," iterations."
    write(*,*)
    write(*,*) "what is your initial guess for the upper limit at this percentile?" 
    read(*,*) ul

    do 
        !Call subroutine with romberg
        call romberg(ll,ul,f,valu,epsi,maxit,n)
        if (abs(valu-percent) < epsi) then
            perc = valu*100
            write(*,'(A,f8.4,A)') "Upper limit reached at ",perc," %"
            max_flow = ul
            exit
        else if (valu>percent) then
            !write(*,'(A,f8.6,A)') "Integration result ",valu," is too high."
            if (abs(valu-percent)>0.05) then
            ul = ul - 4
            else if (abs(valu-percent)<0.05 .or. abs(valu-percent)>0.01) then
            ul = ul - .1
            else if (abs(valu-percent)<0.01) then
            ul = ul - .0001
            end if
        else if (valu<percent) then
            !write(*,'(A,f8.6,A)') "Integration result ",valu," is too low."
            if (abs(valu-percent)>0.05) then
            ul = ul + 3
            else if (abs(valu-percent)<0.05 .or. abs(valu-percent)>0.01) then
            ul = ul + .1
            else if (abs(valu-percent)<0.01) then
            ul = ul + .0001
            end if
        end if
    end do 

    write(*,'(A,f10.6,A)') "Maximum design flow is ", max_flow, " (m^3/day)"
    write(*,'(A,I3)') "Iterations necessary to converge", n

    stop
end program romberg_integration
!________________________________________________________________________

function f(x)
    use vars
    implicit none
    real(wp), intent(in) :: x
    real(wp) :: f, one, two, three

    one = 1/(alpha*x*Gamma(beta))
    two = ((log(x) - lambda)**(beta-1))/(alpha**(beta-1))
    three = exp(-1*((log(x)-lambda)/alpha))
   
    f = one*two*three
    
    return
end function f
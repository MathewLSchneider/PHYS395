! Hw5_1.f90 
! compile with: gfortran -O3 -fdefault-real-8 A5Q2_even.f90 

program Q1
implicit none

real, parameter :: hbar = 1.0, m = 1.0
real psi1, psi2, psi3, psi4, E1, E2, E3
real E_limit1(5), E_limit2(5)

integer, parameter :: iterations = 50

real a, b, c, fa, fb, fc
integer i, j, n
real(8), allocatable :: data(:,:)
real EigenV(5)

n = floor(8.0/1e-3)
! Pre-estimating the bracket for each 5 of the even solutions
do i=1,5
    E_limit1(i)=0.1+2.0*(i-1.0);
    E_limit2(i)=1.0+2.0*(i-1.0);
end do

!E_limit1(1) = 0.3
!E_limit2(1) = 1.0


allocate(data(5,n))

do j=1,5
        E1=E_limit1(j); E2=E_limit2(j);
        E3=bisection(E1,E2, 1.0, (-1.0)**(j-1.0)); !subroutine 
        EigenV(j) = E3
        ! fa=integrate(a, guess, neg*x, neg*dx,j,E3) 
        ! fb=integrate(a, guess, x, dx,(j+1),E3) 
        data(j,:)=integrate((-1.0)**(j-1.0), 8.0, 1e-3, E3, 2.0) 
end do 

write(*,*) "The eigenvalues for the even solutions are...."
do i=1,5 
  write(*,*) EigenV(i)
end do 
!fa = integrate(0.70710678119832193, 10.0, 0.01); stop "Done"

! initial interval
! a = 0.3
! b = 1.0
! fa = integrate(1.0, 8.0, 1e-3, a)
! fb = integrate(1.0, 8.0, 1e-3, b)

! write(*,*) fa, fb
! if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

! ! bisect the interval
! do i = 1,iterations
!   c = (a+b)/2.0; fc = integrate(1.0, 8.0, 1e-3, c)
!   if (fc == 0.0) exit
  
!   ! Ridder's variation on the basic method
!   c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb)
!   fc = integrate(1.0, 8.0, 1e-3, c)
  
!   write (*,*) c, fc

!   if (fa*fc < 0.0) then; b = c; fb = fc; end if
!   if (fc*fb < 0.0) then; a = c; fa = fc; end if
! end do

! write(*,*) "The eigenvalue is", c 


! !psi1 = integrate(1.0, 5.0, 1e-3)
! !psi2 = integrate(0.0, 0.0, 1.0)
! !psi3 = integrate(0.0, -1.0, 0.0)
!psi4 = integrate(0.0, 0.0, -1.0)

contains

! potential
pure function V(x); intent(in) x
	real V, x, k
	k = 1.0

	V = 0.5*k*x*x
end function

! evaluate derivatives
subroutine evalf(u, dudx, E)
        real u(4), dudx(4), E
        
        associate( x => u(1), psi => u(2), Dpsi => u(3), N => u(4))
        dudx(1) = 1.0; dudx(2) = Dpsi
        dudx(3) = -(E - V(x))*psi*2.0*m/hbar
        dudx(4) = abs(u(2))**2 
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt, E)
        integer, parameter :: s = 5, n = 4
        real y(n), g(n,s), dt, E; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape((/ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1/), (/s,s/))
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i), E)
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

!Bisection_function
function bisection(a,b, num, psi); !intent(in) a,b,E
        real a,b,E,c,fa,fb,fc, bisection, num, psi
            fa = f(a, num, psi)
            fb = f(b, num, psi)    
    
            
            if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."
    
            ! bisect the interval
        do i = 1,iterations
            c = (a+b)/2.0; fc = f(c, num, psi); if (fc == 0.0) exit
            
            ! Ridder's variation on the basic method
            c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c, num, psi)
            
            ! write (*,*) c, fc
            
            if (fa*fc < 0.0) then; b = c; fb = fc; end if
            if (fc*fb < 0.0) then; a = c; fa = fc; end if
        end do
        bisection = c
end function bisection

! integrate equations of motion for a given amount of x (written as t)
function integrate(psi, t, dt, E, num)
	real psi, t, dt, integrate, E, num
	real u(4)
  integer i, n
	
	! start from zero at a given gradient
	u = [0.0, psi, 0.0, 1.0]; 
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt, E); if (abs(u(2)) > 10.0) exit
    if (num == 2) then 
		  !write (*,'(5g24.16)') i*dt, u(1), u(2), u(3), u(4)
      write (*,'(5g24.16)') u(1), u(2), u(2)**2/u(4)/2
    end if 
	end do
	
  if (num == 2) then 
      write(*,*) ""
      write(*,*) ""
  end if 

	call gl10(u,t-n*dt, E)
	
	! return state at time t
	integrate = u(2)
end function

function f(E, num, psi); intent(in) E, num
    real f, E, num, psi
    f = integrate(psi, 8.0, 1e-3, E, num) 
end function

end program Q1
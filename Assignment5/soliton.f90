! soliton.f90  -  solve soliton BVP using shooting method
! compile with: gfortran -O3 -fdefault-real-8 soliton.f90 

program soliton
implicit none

! iterations for bisection cycle
integer, parameter :: iterations = 100

real, parameter :: pi = 3.1415926535897932384626433832795028842Q0

real a, b, c, fa, fb, fc, En
integer i

!fa = integrate(0.70710678119832193, 10.0, 0.01); stop "Done"

! initial interval
!a = 0.01; fa = f(a, 0.5)
!b = 1000.35; fb = f(b, 0.5)

!write (*,*) a, fa
!write (*,*) b, fb

!if (fa*fb > 0.0) stop "The root is not bracketed, bailing out..."

! bisect the interval
!do i = 1,iterations
!	c = (a+b)/2.0; fc = f(c, 0.5); if (fc == 0.0) exit
!	
!	! Ridder's variation on the basic method
!	c = c + (c-a)*sign(1.0,fa-fb)*fc/sqrt(fc*fc-fa*fb); fc = f(c, 0.5)
!	
!	write (*,*) c, fc
!	
!	if (fa*fc < 0.0) then; b = c; fb = fc; end if
!	if (fc*fb < 0.0) then; a = c; fa = fc; end if
!end do

fa=f(0.5,0.5)

contains

! potential
pure function V(x); intent(in) x
	real V, x
	
	V = 0.5*x*x 
end function

! potential derivative
pure function DV(psi); intent(in) psi
	real DV, psi
	
	DV = (psi*psi - 1.0)*psi
end function

! total energy of the dynamical system
function E(u, En); 
        real u(4), E, En
        
        associate( psi => u(1), Dpsi => u(2), x => u(3), Norm => u(4) )
        E = Dpsi**2/2.0 - V(psi)
        end associate
end function

! evaluate derivatives
subroutine evalf(u, dudt, En)
        real u(4), dudt(4), En
        
        associate( psi => u(1), Dpsi => u(2), x => u(3), Norm => u(4) )
        dudt(1) = Dpsi; dudt(2) = 2.0*V(x)*psi - 2.0*En*psi; 
        dudt(3) = 1.0 ; dudt(4) = (psi*psi)
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt, En)
        integer, parameter :: s = 5, n = 4
        real y(n), g(n,s), dt, En; integer i, k
        
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
                        call evalf(y + g(:,i)*dt, g(:,i), En)
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of time
function integrate(psi, t, dt, En)
	real psi, t, dt, integrate, En
	real u(4), E0; integer i, n
	
	! start from zero at a given gradient
	u = [psi, 0.0, 0.0, 1.0]; E0 = En
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt, En); if (abs(u(2)) < 0.00001) exit
    write (*,'(4g24.16)') -i*dt, u(1), u(2), u(4)
		write (*,'(4g24.16)') i*dt, u(1), u(2), u(4)
	!write (*,*) u
	end do
	
	call gl10(u,t-n*dt, En)
	!write (*,*) u
	! return state at time t
	integrate = u(1)/u(4)
end function

! function to find a root of...
function f(psi, En); 
	real f, psi, En
	f = integrate(psi, 50.0, 0.001, En)
  !write (*,*) ''; write (*,*) ''
end function


end program
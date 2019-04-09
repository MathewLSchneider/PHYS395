! Hw5_1.f90 
! compile with: gfortran -O3 -fdefault-real-8 A5Q1_1.f90 

program Q1
implicit none

real, parameter :: hbar = 1.0, m = 1.0, E = 0.0
real psi1, psi2, psi3, psi4

psi1 = integrate(1.0, 5.0, 1e-3)
!psi2 = integrate(0.0, 0.0, 1.0)
!psi3 = integrate(0.0, -1.0, 0.0)
!psi4 = integrate(0.0, 0.0, -1.0)

contains

! potential
pure function V(x); intent(in) x
	real V, x, k
	k = 1.0

	V = 0.5*k*x*x
end function

! potential derivative
pure function DV(x); intent(in) x
	real DV, x, k
	k=1.0

	DV = k*x
end function

! evaluate derivatives
subroutine evalf(u, dudx)
        real u(3), dudx(3)
        
        associate( x => u(1), psi => u(2), Dpsi => u(3))
        dudx(1) = 1.0; dudx(2) = Dpsi
        dudx(3) = -(E - V(x))*psi*2.0*m/hbar 
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dt; integer i, k
        
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
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

! integrate equations of motion for a given amount of x (written as t)
function integrate(pi, t, dt)
	real pi, t, dt, integrate
	real u(3), E0; integer i, n
	
	! start from zero at a given gradient
	u = [0.0, 1.0, 0.0]; E0 = E
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt); if (abs(u(2)) > 10.0) exit
		write (*,'(4g24.16)') i*dt, u(1), u(2), u(3)
	end do
	
	call gl10(u,t-n*dt)
	
	! return state at time t
	integrate = u(1)
end function

end program
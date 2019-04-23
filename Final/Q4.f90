! A4.f90  -  Double Pendulum
! compile with: gfortran -O3 -fdefault-real-8 -fopenmp Q4.f90 -lcfitsio

program Q4
implicit none

! coupling parameter in the potential
real, parameter :: Ei = -1.0 !Here l = g and m = 2
real, parameter :: pi = 3.14159265358979323846264338327950288

! image size and scan bounds
integer, parameter :: nx = 3 * 2**5, ny = 2**6
real(8), parameter :: xx(2) = [-1.5, 1.5], yy(2) = [-1.0, 1.0]

! data array (allocatable to avoid problems with system resources)
real(4), allocatable :: data(:,:,:)
real traj(2)
! temporary variables
integer i, j 
real dt, init

allocate(data(4,nx,ny))

! timestep resolving fastest timescale in the scan
dt = 0.0001

do i = 1,20
    traj = integrate(i/10.0,0.2, 20.0, dt)
end do
do i = 21,29
    traj = integrate(i/10.0,0.3, 20.0, dt)
end do
dt = 0.01
do i = 30,40
    traj = integrate(i/10.0,1.0, 20.0, dt)
end do
dt = 0.1
do i = 41,50
    traj = integrate(i/10.0,5.0, 100.0, dt)
end do
! write out image to file
!call write2fits('data.fit', data, xx, yy, ['x','y','v','w'], '(x0,y0)')

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!dynamical system and Gauss-Legendre integrator!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! potential for 2D particle motion
pure function V(x); intent(in) x
	real V, x
	
	V = -12.0/(cosh(x)*cosh(x))
end function

! total energy of the dynamical system
pure function E(u); intent(in) u
        real u(2), E

        E =  0.5*u(2)*u(2) + V(u(1))
end function

! evaluate derivatives
subroutine evalf(u, dudt)
        real u(2), dudt(2) 


        associate( x => u(1), v => u(2))
        dudt(1) = u(2)		
        dudt(2) = -24.0*sinh(x)/(cosh(x)**3)		

        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 2
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

! integrate equations of motion for a given amount of time
function integrate(x, init, t, dt)
	real x, t, dt, integrate(2)
	real u(2), E0, init; integer i, n

	! start from a given position at rest
	u = [x, 0.0]; E0 = V(x)
	
	! number of time steps needed
	n = floor(t/dt)
	
	do i = 1,n
		call gl10(u, dt); 
        if (i*dt > init) then
            if (abs(u(2)) < 0.01 ) exit;
        end if
		!write (*,'(4g24.16)') i*dt, u(1), u(2), E(u)/E0 - 1.0
		                            
	end do
	
    write (*,'(4g24.16)') 2.0*i*dt, E(u)
	
	! return state at time t
	integrate = u
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write array data into FITS file as sequence of image extensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write2fits(file, array, xx, yy, vars, coords)
        character(len=*) file, vars(:), coords
        real(4) array(:,:,:); real(8) xx(2), yy(2)
        optional xx, yy, vars, coords
        
        integer i, j, status, unit
        integer :: hdus, naxis = 2, n(2), npix
        integer :: bitpix = -32, group = 1, blocksize = -1
        
        ! data dimansions
        hdus = size(array,1)
        n(1) = size(array,2)
        n(2) = size(array,3)
        npix = n(1)*n(2)
        
        ! delete file if it already exists
        open(unit=1234, iostat=status, file=file, status='old')
        if (status == 0) close(1234, status='delete'); status = 0
        
        ! initialize FITS file
        call ftgiou(unit, status)
        call ftinit(unit, file, blocksize, status)
        
        ! write image extensions
        do i = 1,hdus
                call ftiimg(unit, bitpix, naxis, n, status)
                call ftppre(unit, group, 1, npix, array(i,:,:), status)
                
                if (present(vars)) then
                        if (present(coords)) then
                                call ftpkys(unit, 'EXTNAME', trim(vars(i))//coords, 'variable stored in extension', status)
                        else
                                call ftpkys(unit, 'EXTNAME', trim(vars(i)), 'variable stored in extension', status)
                        end if
                end if
                if (present(xx)) then
                        call ftpkyj(unit, 'CRPIX1', 1, 'x-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL1', xx(1), 14, 'x-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT1', (xx(2)-xx(1))/(n(1)-1), 14, 'x-axis increment', status)
                end if
                if (present(yy)) then
                        call ftpkyj(unit, 'CRPIX2', 1, 'y-axis origin pixel', status)
                        call ftpkyd(unit, 'CRVAL2', yy(1), 14, 'y-axis origin coordinate', status)
                        call ftpkyd(unit, 'CDELT2', (yy(2)-yy(1))/(n(2)-1), 14, 'y-axis increment', status)
                end if
        end do
        
        ! clean up
        call ftclos(unit, status)
        call ftfiou(unit, status)
end subroutine write2fits

end program

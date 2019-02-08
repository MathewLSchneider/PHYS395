

program leastsq2
implicit none

integer, parameter :: n = 4
real x, y, A(n,n), B(n), U(n), Chi
real, dimension(:), allocatable :: fit, xp

integer i, j
integer :: s = 0
integer :: status = 0
integer :: stat = 0
Chi = 0.0

open (unit=1,file='Fit.txt',action='read', &
	status='old',position='append')
do while (status < 1)
	read (1,*,iostat = status) s, U
	if (status < 0) then
		write (*,*) 'problem'
		stop
	end if  
end do

write (*,*) s
close (1)
open (unit = 2, file='Assignment #2.dat', action='read', &
	status='old', position='rewind')
allocate(fit(s)); allocate(xp(s))
fit = 0; xp = 0
do i=1,s
	read (2,*) x, y
	do j=1,n
		fit(i) = fit(i) + U(j)*cos(2.0*(j-1)*355.0*x/113.0)
	end do 
	xp(i) = x
	Chi = Chi + (fit(i)-y)**2
end do
close(2)

open (unit=3,file='DATA.txt',action='write', &
	status='new',position='append')
do i=1,s-3
	write (1,*) xp(i), fit(i)
end do
close(3)
end program
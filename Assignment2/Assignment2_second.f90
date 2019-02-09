

program Assignment2_second
implicit none

integer, parameter :: n = 5
real x, y, U(n), Chi
real, dimension(:), allocatable :: fit, xp

integer i, j
integer :: s = 0
integer :: status = 0
integer :: stat = 0
Chi = 0.0

!Read in coefficients
open (unit=1,file='Fit.txt',action='read', &
	status='old')
do i=1,1
	read (1,*) U
end do
s = U(1)
close (1)

!Read in data
open (unit = 2, file='Assignment #2.dat', action='read', &
	status='old')
allocate(fit(s)); allocate(xp(s))
fit = 0; xp = 0
do i=1,s
	read (2,*) x, y
	do j=1,n-1
		fit(i) = fit(i) + U(j+1)*cos(2.0*(j-1)*355.0*x/113.0)
	end do 
	xp(i) = x
	Chi = Chi + (fit(i)-y)**2
end do
close(2)

!Evaluate goodness of fit
write (*,*) 'Chi squared is ', Chi
if (Chi < s + s/10 .AND. Chi > s - s/10) then
	write (*,*) 'Chi squared is comparable to the number of points',&
			s, '.  Fit looks good'
else if (Chi < s - s/10) then
	write (*,*) 'Chi squared is a bit low, perhaps overfitted'
else if (Chi > s + s/10) then
	write (*,*) 'Chi squared is a bit high, maybe refit'
end if

!Write fit to data file
open (unit=3,file='DATA.txt',action='write', &
	status='new',position='append')
do i=1,s
	write (3,*) xp(i), fit(i)
end do
close(3)


end program
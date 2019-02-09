

program Assignment2_fifth
implicit none

integer, parameter :: n = 4
real x, y, U(n), Chi, B(n), Sing(n)
real, dimension(:), allocatable :: fit, xp, W, D

real, dimension(:,:), allocatable :: A, F

integer i, j, r, info
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

!Read in the data
open (unit = 2, file='Assignment #2.dat', action='read', &
	status='old',position='rewind')
allocate(A(s,n)); allocate(xp(s)); allocate(F(s,1)); allocate(D(s))
A = 0; F = 0; B = 0
do i=1,s
	read (2,*) x, y
	call evalb(x,B)
	A(i,:) =  B !Form the matrix for lss
	xp(i) = x
	F(i,1) = y 
	D(i) = y
end do
close(2)

allocate(W(s))

!Solve the lss problem with dgelss
call dgelss(s,3,1,A,s,F,s,Sing,-1.0,r,W,s,info)

write (*,*) 'The condition number is:'
write (*,*) maxval(Sing)/minval(Sing, mask = Sing .GT.0)

!Get the coefficients
U(:) = F(1:n,1)
write (*,*) 'Coefficients using LAPACK dgelss are:'
write (*,*) U

!Create the fit line
allocate(fit(s)); fit = 0
do i=1,s
	do j=1,n
		fit(i) = fit(i) + U(j)*cos(2.0*(j-1)*xp(i)*355.0/113.0)
	end do 
end do		


Chi = sum((fit-D)**2)

!Check goodness of fit
write (*,*) 'Chi squared is ', Chi
if (Chi < s + s/10 .AND. Chi > s - s/10) then
	write (*,*) 'Chi squared is comparable to the number of points',&
			s, '.  Fit looks good'
else if (Chi < s - s/10) then
	write (*,*) 'Chi squared is a bit low, perhaps overfitted'
else if (Chi > s + s/10) then
	write (*,*) 'Chi squared is a bit high, maybe refit'
end if

!Write fit to file
open (unit=3,file='DATA.txt',action='write', &
	status='old',position='rewind')
do i=1,s
	write (3,*) xp(i), fit(i)
end do
close(3)

contains

subroutine evalb(x, B)
	real x, B(n)
		
	forall (i=0:n-1) B(i+1) = cos(2.0*i*x*355.0/113.0)
end subroutine

end program
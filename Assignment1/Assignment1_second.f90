!100 term approx.

program hundredterms
implicit none

integer  n, i, j
real F(100), G(100), x(100), T(100,100), U(100,100), C(100), D(100)


n = 100

!Initializing
call Initialize(n,x,T,U,F,G)

!Solve for coeff.
call gaussj(n, T, F)

!Sum over basis
do i=1,n
	do j = 1,n
		C(i) = C(i) + F(j)*ChebyshevT(x(i),j-1)
		D(i) = D(i) + F(j)*DerChebyshevT(x(i),j-1)
	end do
end do

!Black magic, do not touch
C(79) = 0
C(71) = 0
do j = 1,n
	C(79) = C(79) + F(j)*ChebyshevT(x(79), j-1)
end do
do j = 1,n
	C(71) = C(71) + F(j)*ChebyshevT(x(71), j-1)
end do

!Print function
do i=1,n
	write (*,*) x(i), (1.0/(1.0 + 10.0*x(i)*x(i)))
end do

write (*,*) ""
write (*,*) ""

!Print approx.
do i=1,n
	write (*,*) x(i),  C(i)
end do

write (*,*) ""
write (*,*) ""

forall (i=1:n) F(i) = (-20.0*x(i)/(1.0 + 10.0*x(i)*x(i))**2)

!Print deriv.
do i=1,n
	write (*,*) x(i), F(i)
end do 

write (*,*) ""
write (*,*) ""

!Print deriv.
do i=1,n
	write (*,*) x(i),  D(i)
end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!Initialize matrices
subroutine Initialize(n,x,T,U,F,G)
	integer n, i, j
	real x(n), T(n,n), F(n), U(n,n), G(n)

	do i = 1,n
		x(i) = -1.0 + (i-1.0)*2.0/(n-1)
	end do

	forall (i=1:n) F(i) = (1.0/(1.0 + 10.0*x(i)*x(i)))
	forall (i=1:n) G(i) = (1.0/(1.0 + 10.0*x(i)*x(i)))

	do i=1,n
		do j=1,n
			T(i,j) = ChebyshevT(x(i),j-1)
		end do
	end do

	do i=1,n
		do j=1,n
			U(i,j) = DerChebyshevT(x(i),j)
		end do
	end do

end subroutine

! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n, i, j; real A(n,n), B(n), M(n,n+1)

!Make matrix M which is just [A|B]
	M(:,1:n) = A(:,:)
	M(:,n+1) = B(:)

	!Triangularize the A part of the matrix
	do j = 1,n-1
		do i = j+1,n
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
	end do

	!Diagonalize and normalize the A part of the matrix
	do j = n,2,-1
		do i = j-1,1,-1
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
		M(j,:) = M(j,:)/M(j,j)
	end do

	
	!The last column is now the solution
	B(:) = M(:,n+1)

end subroutine


! Chebyshev polynomial T_n(x)
elemental function ChebyshevT(x, n)
	real ChebyshevT, x; integer n
	intent(in) x, n
	
	ChebyshevT = cos(n*acos(x))
end function

!Deriv. of Cheb.
elemental function DerChebyshevT(x, n)
	real DerChebyshevT, x; integer n
	intent(in) x, n
	
	DerChebyshevT = n*sin(n*acos(x))/sqrt(1.0-x*x)
end function

end
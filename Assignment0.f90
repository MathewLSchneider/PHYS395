! compile with: gfortran Assignment0.f90 -o Fibonacci
! run with:     ./Fibonacci > DATA


program Fibonacci
implicit none

integer*8 a, b, c, i !make the variables 8-bit integers as lots of digits are needed
integer, parameter :: n = 48 !counting parameter, 2 off since we start with F0 and F1

!Initialize first two values of Fibonacci
a = 0
b = 1

write (*,*) a
write (*,*) b !Print the two first numbers

do i = 1,n

  c = a + b !!Generate the next Fibonacci number
  a = b     !!Shift over the variable a and b
  b = c
  write (*,*) c  !Print only the next number

end do

end

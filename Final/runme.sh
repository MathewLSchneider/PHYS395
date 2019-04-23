#!/bin/bash
chmod u+x plot0.gpl plot1.gpl plot3.gpl

echo "Doing question 1 and 2..."
gfortran -O3 -fdefault-real-8  Q1_2.f90 
./a.out

echo "Doing question 3"
gfortran -O3 -fdefault-real-8 -fopenmp Q3.f90 -lcfitsio
./a.out > DATA

./plot0.gpl

echo "Doing question 4"
gfortran -O3 -fdefault-real-8 -fopenmp Q4.f90 -lcfitsio
./a.out > DATA

./plot1.gpl

echo "Doing question 5"
gfortran -O3 -fdefault-real-8 Q5first.f90 -llapack
./a.out

gfortran -O3 -fdefault-real-8 Q5second.f90 -llapack
./a.out > DATAeven

gfortran -O3 -fdefault-real-8 Q5third.f90 -llapack
./a.out > DATAodd

./plot3.gpl
echo "Plots show the even then odd solutions"
echo "Energy eigenvalues are in EigVal.txt"
echo "All Done"
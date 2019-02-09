#!/bin/bash
#to run, first enter chmod u+x runme.sh plot1.gpl plot2.gpl plot3.gpl
#then enter ./runme.sh
#If you want to run again, first delete the Fit.txt and DATA.txt files

echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo 'Starting with n=3:'
gfortran -O3 -fdefault-real-8 -o A2 Assignment2_first.f90 -llapack
./A2 < 'Assignment #2.dat'
echo "Fitting done for n=3 coefficients"
gfortran -O3 -fdefault-real-8 -o A2 Assignment2_second.f90 -llapack
./A2 


./plot1.gpl

echo 'n = 3 plotted'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

echo 'Now doing n=8'
gfortran -O3 -fdefault-real-8 -o A2 Assignment2_third.f90 -llapack
./A2 < 'Assignment #2.dat'
echo "Fitting done for n=8 coefficients"
gfortran -O3 -fdefault-real-8 -o A2 Assignment2_fourth.f90 -llapack
./A2 

./plot2.gpl

echo 'n=8 plotted'
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo 'Now using LAPACK...'
gfortran -O3 -fdefault-real-8 -o A2 Assignment2_fifth.f90 -llapack
./A2 
echo "Fitting done using LAPACK dgelss"

./plot3.gpl


echo "All Done"

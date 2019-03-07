#!/bin/bash
#to run, first enter chmod u+x runme.sh plot1.gpl 
#then enter ./runme.sh
#If you want to run again, first delete the Fit.txt and DATA.txt files

echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo "Doing questions 1 through 3..."
gfortran -O3 -fdefault-real-8 -o A3 A3Q1_3.f90
./A3
echo "Now moving to 4 and 5..."
gfortran -O3 -fdefault-real-8 -o A3 A3Q4_5.f90 -llapack
echo "Fitting data..."
./A3 < 'Assignment #2.dat'

echo "Plotting data and fit..."
./plot1.gpl

echo "Fit saved"

echo "All Done"

#!/bin/bash
#to run, first enter the following:
#chmod u+x runme.sh animating.gpl save.py savezoom1.py savezoom2.py savezoom3.py savezoom4.py
#then enter ./runme.sh
chmod u+x plot0.gpl plot1.gpl plot2.gpl plot3.gpl plot4.gpl 

echo "Showing bad eigenvalue"
gfortran -O3 -fdefault-real-8 -o A5 A5Q1_1.f90 
./A5 > DATAeven

gfortran -O3 -fdefault-real-8 -o A5 A5Q1_2.f90 
./A5 > DATAodd

./plot0.gpl

echo "Shooting for QSHO solution"
gfortran -O3 -fdefault-real-8 -o A5 A5Q2_even.f90 
./A5 > DATAeven


gfortran -O3 -fdefault-real-8 -o A5 A5Q2_odd.f90 
./A5 > DATAodd

echo "Plotting Wavefunction"
./plot1.gpl

echo "Shooting for Anharmonic Solution"
gfortran -O3 -fdefault-real-8 -o A5 A5Q3_even.f90 
./A5 > DATAeven

gfortran -O3 -fdefault-real-8 -o A5 A5Q3_odd.f90 
./A5 > DATAodd

echo "Plotting Wavefunction"
./plot2.gpl

echo "Relaxation for QSHO"
gfortran -O3 -fdefault-real-8 rayleigh4.f90 -llapack
./a.out

gfortran -O3 -fdefault-real-8 integrating1.f90 
./a.out > DATAeven
gfortran -O3 -fdefault-real-8 integrating_odd1.f90 
./a.out > DATAodd

echo "Plotting"
./plot3.gpl

gfortran -O3 -fdefault-real-8 rayleigh5.f90 -llapack
./a.out

echo "Relaxtion for anharmonic"
gfortran -O3 -fdefault-real-8 integrating.f90 
./a.out > DATAeven
gfortran -O3 -fdefault-real-8 integrating_odd.f90 
./a.out > DATAodd

echo "Plotting"
./plot4.gpl

echo "All Done"
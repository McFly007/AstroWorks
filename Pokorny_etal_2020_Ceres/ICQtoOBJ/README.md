Written by Petr Pokorny on November 23, 2020
----
This is a simple FORTRAN code inspired by https://sbnarchive.psi.87edu/pds3/dawn/fc/DWNCSPC401/DOCUMENT/ICQMODEL.ASC
for conversion of the ICQ format to the OBJ format. 

This folder contains two files:
(1) 01_runme.sh - a bash script that downloads the ICQ data file for Ceres used in our paper, compiles, and runs the code
(2) ICQtoOBJ.f90 - the FORTRAN code
(2) OBJ2BINARY.f90 - the FORTRAN code

What to do:
Open the terminal/command line and write: bash 01_runme.sh or make 01_runme.sh executable and run it.
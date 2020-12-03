#!/bin/bash

# Get the ICQ topography file from the DAWN archive
wget 'https://sbnarchive.psi.edu/pds3/dawn/fc/DWNCSPC_4_01/DATA/ICQ/CERES_SPC181019_0512.ICQ'

# Compile the code for ICQ to OBJ conversion
gfortran ICQtoOBJ.f90 -o ICQtoOBJ -O3

gfortran OBJ2BINARY.f90 -o OBJ2BINARY -O3

# Runs the code - the outputs are 
# (1) CERES_0512.obj the obj data file
# (2) CERES_0512.info that containts number of vectors and faces of the obj data file - you can check the number of vectors and faces it by running awk for example: awk '{if($1=="f") {f+=1}; if($1=="v") {v+=1}} END {print f,v}' CERES_0512.obj
./ICQtoOBJ


# Converts the OBJ file to a FORTRAN UNFORMATTED file that is faster to read
./OBJ2BINARY



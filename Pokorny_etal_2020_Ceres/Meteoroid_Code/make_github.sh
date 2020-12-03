#!/bin/sh

#CXXFLAGS := -Wall -Wextra -Werror -Wfatal-errors
#CXXFLAGS := $(CXXFLAGS) -I $(CURDIR)/include 
#CXXFLAGS := $(CXXFLAGS) -std=c++14
CDIR=./C_codes

g++ $CDIR/BBox.cpp -o $CDIR/BBox.o -c  -Wall -Wextra -msse3 -Ofast
g++ $CDIR/BVH.cpp -o $CDIR/BVH.o -c -Wall -Wextra -msse3 -Ofast
g++ $CDIR/simple_vis.cpp -o $CDIR/simple_vis.o  -c -Wall -Wextra -msse3 -Ofast
gfortran RAY_MODS_Github.f90 -o RAY_MODS_Github.o -c -Ofast -fcheck=all -fbacktrace
gfortran $CDIR/simple_vis.o $CDIR/BBox.o $CDIR/BVH.o RAY_MODS_Github.o Meteoroids_CODE_V03_Github.f90 -lstdc++ -o Meteoroids_CODE_V03_Github -Ofast -fopenmp -fcheck=all -fbacktrace

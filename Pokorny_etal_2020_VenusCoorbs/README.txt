Files for Pokorny, Kuchner and Sheppard (2020)
----
These are the files for Pokorny, Kuchner and Sheppard (2020). Each file in this directory is described below. This directory will be updated after the paper is reviewed by the referee, so we can address any requested changes. The code "Asteroid_Types.f90" is a Fortran code that calculates the visual magnitudes of asteroids as observed from Earth using Penttila et al. (2016) formalism and G1 and G2 values from Shevchenko et al. (2016).
----
(1) Asteroid_Types.f90 - The main code that reads file "Datafile.dat" and outputs Visual_Magnitudes
(2) Datafile.dat - contains a sample asteroid population - column description (1 - time, 2 - particle id, 3- Heliocentric Distance (au), 4- Distance From Earth (au), 5 - Phase Angle (radians), 6 - Heliocentric Ecliptic Longitude (radians), 7 - Heliocentric Ecliptic Latitude (radians)) - only #3, #4, #5 are necessary for the calculation
(3) Visual_Magnitudes - contains visual magnitudes for D = 0.5, 1.0, 1.5, 2.0, 2.5 km and 6 asteroid types S, M, E, C, P and D = 30 columns

Read me for Meteoroid Bombardment of Lunar Poles by Pokorny et al. (2020)

----- COMPILING/RUNNING THE PROGRAM -----
(1) Run ./@make to compile impact_ray_simul_V002_GITHUB.f90 - the @make script is primitive, essentially just calls "gfortran impact_ray_simul_V002_GITHUB.f90 -o impact_ray_simul_V002_GITHUB -Ofast -mcmodel=medium"
(2a) Run ./impact_ray_simul_V002_GITHUB < C_NORTH.in - for the North Pole calculation 
(2b) Run ./impact_ray_simul_V002_GITHUB < C_SOUTH.in - for the South Pole calculation
(3) The program writes a binary file containing 6 values per line : x,y,z of the first vectex of a triangle face (3 values), meteoroid mass flux in (g/cm^2/s), meteoroid energy flux in (kJ/cm^2/s), and meteoroid ejecta production rate (g/cm^2/s) - all values are real*4 (float)
------------------------------------------

----- FILES IN THIS DIRECTORY -----
(1) @make - a simple compiling script in csh
(2) AST - radiant/velocity distribution of main-belt asteroids at Moon averaged over one year. (#1) is the heliocentric ecliptic longitude (deg), (#2) is the heliocentric ecliptic latitude (deg), (#3) is the impact velocity (km/s), (#4) is the meteoroid mass flux in (kg) normalized to 1000 kg per day accretion rate at Earth (i.e. assuming Earth's cross section) - to get units kg/m^2, divide (#4) by 6378^2*1000^2*pi
(3) C_NORTH.in - input file for North Pole calculation - the minimum latitude is set to -15 degrees to shorten the calculation. It is possible to split this calculation into many steps and then add the results afterwards - The entire calculation takes approximately 100 CPU hours
(4) C_SOUTH.in - the same as (3) but for the South Pole
(5) Data_North_Pole.tar.gz - tarball containing an ASCII file with North Pole Results
(6) Data_South_Pole.tar.gz - tarball containing an ASCII file with South Pole Results
(7) GRID_North_Pole.bin - binary gridded Lunar North Pole topography used in our paper derived from LOLA dataset (set the paper for details)
(8) GRID_South_Pole.bin - the same as (7) but for the South Poles
(9) HTC - the same as (2) but for the Halley-type Comets
(10) JFC - the same as (2) but for the Jupiter-family Comets
(11) OCC - the same as (2) but for the Oort Cloud Comets
(12) README.txt - this file
(13) header - describes all columns in (5) and (6)
(14) header_1line - the one line header for (5) and (6)
(15) impact_ray_simul_V002_GITHUB.f90 - the latest version of our code
(16) vars.in - variables definitions for (15)
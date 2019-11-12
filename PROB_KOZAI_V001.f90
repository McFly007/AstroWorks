      PROGRAM KOZAIOPIKPROB
! CHANGE LOG:
! VER 001 - Mar 1 2017 - by Petr Pokorny (petr.pokorny@nasa.gov)
! -----------------------
! o Now everything works for any apl and epl (previously it required some tuning, my bad)
! o For epl < 1e-5 the code checks only impacts at apl while for epl >= 1e-5 the code scans the whole range from perihelion to aphelion of the target
! o The output was redesigned to write for each impact configuration the heliocentric distance of the target, the radiant location (phi, theta), the impact velocity and the impact probability
! o Added the output file header
! o Individual impact probabilities are now properly weighted and their sum gives the total intrinsic impact probability (previously it was not like that)
! o The program writes the orbital elements of the projectile and its total impact probability with the target on the screen
! -----------------------
! Author:
! Petr Pokorny
! 
! LICENCE:
!                 ____________  
! This program is !!! FREE !!! to use, however, it would be very kind of you to
!                 ------------
! put a reference into your article (reference to Vokrouhlicky et al. 2012)
! or better an article with Pokorny, Vokrouhlicky 2013, where is the code cited (^_^)


! Reference publication: 
!
! Opik-type collision probability for high-inclination orbits
! Vokrouhlicky, David; Pokorny, Petr; Nesvorny, David (2012)
! Icarus, Volume 219, Issue 1, p. 150-160.
!
! Opik-type collision probability for high-inclination orbits: Targets on eccentric orbits
! Pokorny, Petr; Vokrouhlicky, David (2013)
! Icarus, Volume 226, Issue 1, p. 682-693

! SUBROUTINE KOZAIOPIK computes impact probabilities between Target and Projectile
! The Target is assumed to be in eccentric and coplanar orbit (i = 0, 0 =< e < 1)
! The Projectile is assumed to be on bound elliptic orbit around the Sun and non-resonant with any other planet
! The Projectile is defined by its semimajor axis (AU), eccentricity, Inclination (deg) and argument of pericenter (deg)
!
!                                /
!                               /
!                              /                       ___ 
!                             *                    ,o88888   
!                                               ,o8888888'  
!                         ,:o:o:oooo.        ,8O88Pd8888"    
!                     ,.::.::o:ooooOoOoO. ,oO8O8Pd888'"     
!                   ,.:.::o:ooOoOoOO8O8OOo.8OOPd8O8O"     
!                  , ..:.::o:ooOoOOOO8OOOOo.FdO8O8" 
!                 , ..:.::o:ooOoOO8O888O8O,COCOO" 
!                , . ..:.::o:ooOoOOOO8OOOOCOCO" 
!                 . ..:.::o:ooOoOoOO8O8OCCCC"o 
!                    . ..:.::o:ooooOoCoCCC"o:o 
!                    . ..:.::o:o:,cooooCo"oo:o: 
!                 `   . . ..:.:cocoooo"'o:o:::' 
!                 .`   . ..::ccccoc"'o:o:o:::' 
!                :.:.    ,c:cccc"':.:.:.:.:.' 
!              ..:.:"'`::::c:"'..:.:.:.:.:.' 
!            ...:.'.:.::::"'    . . . . .' 
!           .. . ....:."' `   .  . . '' 
!         . . . ...."' 
!         .. . ."'    
!        . 
!
!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      INTEGER nmax ! Maximal number of iterations in the evaluation of the integral
      INTEGER coef ! Maximal number of steps used for the evaluation in all iterations - defines output arrays
      parameter(nmax=7) ! Determines the precision of the evaluation of the integral. The value 7 seems to give the best ratio between speed and precision
      parameter(coef=2*2*(3**(nmax-1))) ! Formula for the total number of steps of the integral

      real*8 outP(coef,8) ! All intrinsic collisional probabilities [au^-2 yr^-1] evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outTH(coef,8) ! All theta [deg] of the radiants evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outPH(coef,8) ! All phi [deg] of the radiants evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outR(coef) ! All heliocentric distances [AU] of the target evaluated in the integral
      real*8 outVEL(coef,8) ! All impact velocities evaluated in the integral | 8 possible impact configurations for one r_0


      real*8 apl     !   Semimajor Axis of the target [AU]
      real*8 epl     !   Eccentricity of the target

      real*8 a       !   Semimajor Axis of the projectile [AU]
      real*8 e       !   Eccentricity of the projectile
      real*8 inc     !   Inclination of the projectile [deg]
      real*8 omega   !   Argument of pericenter of the projectile [deg]
      real*8 outres  !   Total intrinsic collisional probability in au^-2 yr^-1
      real*8 vel(8)
      real*8 lambda(8)
      real*8 beta(8)
      real*8 prob(8)
      integer j,l
      integer cntr
      
      real*8 pi       
      data pi/3.141592653589793/


!     write(*,*) "Semimajor axis of projectile in [AU]"
!     read(*,*)
!     do
      read(*,*) a,e,inc,omega

! !     write(*,*) "Eccentricity of projectile"
!       read(*,*) e
! 
! !     write(*,*) "Inclination of projectile in (deg)"
!       read(*,*) inc
! 
! !     write(*,*) "Argument of projectile pericenter (deg)"
!       read(*,*) omega

!     Degrees to radians
      inc=inc/180.d0*pi
      omega=omega/180.d0*pi


!     write(*,*) "Semimajor Axis of the target [AU]"
      read(*,*) apl
!      apl=1.0
!       epl=0.001671123d0
!      epl=0.016d0

!     write(*,*) "Eccentricity of the target"
       read(*,*) epl

      open(21,file="output.txt", status="unknown")
        write(21,*) "  Target's heliocentric | Longitude of the impact | &
                      Latitude of the impact |  Impact velocity (km/s) |   & 
                        Impact probability    "
        write(21,*) "     distance (au)      |    vs target's apex     |   &
                        vs target's apex     |                         |     & 
                           (au^-2 yr^-1)"
        write(21,*) "--------------------------------------------------------------&
                     --------------------------------------------------------------"
!     Calling the main subroutine
      call WETKOZ(apl,epl,a,e,inc,omega,outres,vel,outP,outTH,outPH,outR,outVEL)
!     This recalculates the collisonal probability for
!     Earth-sized object
!     outres=outres*(6378d0/1.5E8)**2

!     Write an output: Total intrinsic collisional probability
      write(*,*) a,e,inc*180d0/pi,omega*180d0/pi,outres
!       write(*,*) "----------"
      do l=1,coef
!        cntr=0
       do j=1,8
       if (outp(l,j).gt.0) then
        write(21,1000) outR(l),outPH(l,j),outTH(l,j),outVEL(l,j),outp(l,j)
!        cntr=cntr+1
       endif
!        write(*,*)j, prob(j)
       enddo
!       write(*,"(A,F7.3,I1)") "Number of impact configurations for R = ",outR(l),cntr
      enddo
!         write(*,*) lambda
!         write(*,*) beta
!         write(*,*) vel
!         write(*,*) prob
!       write(*,*) "----------"
      close(21)
 1000 FORMAT(3X,5(ES19.8,6X))
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         .--..--..--..--..--..--.
!       .' \  (`._   (_)     _   \
!     .'    |  '._)         (_)  |
!     \ _.')\      .----..---.   /   SUBROUTINE AND FUNCTIONS AREA
!     |(_.'  |    /    .-\-.  \  |
!     \     0|    |   ( O| O) | o|
!      |  _  |  .--.____.'._.-.  |
!      \ (_) | o         -` .-`  |
!       |    \   |`-._ _ _ _ _\ /
!       \    |   |  `. |_||_|   |
!       | o  |    \_      \     |     -.   .-.
!       |.-.  \     `--..-'   O |     `.`-' .'
!     _.'  .' |     `-.-'      /-.__   ' .-'
!   .' `-.` '.|='=.='=.='=.='=|._/_ `-'.'
!   `-._  `.  |________/\_____|    `-.'
!      .'   ).| '=' '='\/ '=' |
!      `._.`  '---------------'
!              //___\   //___\
!                ||       ||
!                ||_.-.   ||_.-.
!               (_.--__) (_.--__)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      SUBROUTINE WETKOZ(apl,epl,apr,e,inc,omega,outres,outv,outP,outTH,outPH,outR,outVEL)

      IMPLICIT NONE

      INTEGER n ! Counter for the number of iterations in computing of the integral
      INTEGER it ! Number of steps in the integral
      INTEGER j ! Loop integer
      INTEGER nmax ! Maximal number of iterations in the evaluation of the integral
      INTEGER coef ! Maximal number of steps used for the evaluation in all iterations - defines output arrays
      INTEGER l ! Counter determining the position of the element in arrays
      INTEGER k ! Loop integer
      INTEGER ik ! Integer 1 or 2 determining whether the true anomaly of the target is (0,180) or (180,360) - we to investigate the whole orbit

      parameter(nmax=7) ! Determines the precision of the evaluation of the integral. The value 7 seems to give the best ratio between speed and precision
      parameter(coef=2*2*(3**(nmax-1))) ! Formula for the total number of steps of the integral


!===== INPUTS =====
      real*8 apl    !   Semimajor Axis of the target [AU]
      real*8 epl    !   Eccentricity of the target

      real*8 apr    !   Semimajor Axis of the projectile [AU]
      real*8 e      !   Eccentricity of the projectile
      real*8 inc    !   Inclination of the projectile [radians]
      real*8 omega  !   Argument of pericenter of the projectile [radians]
!===== END OF INPUTS =====


!===== OUTPUTS ======
      real*8 outres ! Output: Total intrinsic collisional probability in au^-2 yr^-1
!===== END OF OUTPUTS =====


      real*8 prob(8)!   Intrinsic collisional probability for 8 possible impact configurations for one r_0 of the target
      real*8 prob1  !   Sum of prob(8), total intrinsic collisional probability for one r_0 of the target


      real*8 pi
      real*8 theta(8) ! Theta of the radiants for 8 possible impact configurations
      real*8 phi(8) ! Phi of the radiants for 8 possible impact configurations
      real*8 outk(8) ! k = e*cos(omega) for 8 possible impact configurations
      real*8 outh(8) ! h = e*sin(omega) for 8 possible impact configurations
      real*8 outv(8) ! Impact velocity



!===== POSSIBLE OUTPUTS =====
      real*8 outR(coef) ! All heliocentric distances [AU] of the target evaluated in the integral
      real*8 outP(coef,8) ! All intrinsic collisional probabilities [au^-2 yr^-1] evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outVEL(coef,8) ! All impact velocities evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outTH(coef,8) ! All theta [deg] of the radiants evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outPH(coef,8) ! All phi [deg] of the radiants evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outKK(coef,8) ! All  k = e*cos(omega) evaluated in the integral | 8 possible impact configurations for one r_0
      real*8 outhh(coef,8) ! All  h = e*sin(omega) evaluated in the integral | 8 possible impact configurations for one r_0
!===== END POSSIBLE OUTPUTS =====



      real*8 a ! Lower limit of the integral
      real*8 b ! Upper limit of the integral
      real*8 x ! r_0 in the integral
      real*8 del ! Integration step (see Press et al. 1992, chap. 4.1)
      real*8 ddel  ! 2 * del
      real*8 s1 ! Intermediate sum of the first integral [a ( 1 - e ) , a ]
      real*8 s2 ! Intermediate sum of the second integral [ a , a ( 1 + e ) ]

      real*8 aa ! Pericenter of the target
      real*8 bb ! Apocenter of the target
      real*8 tnm ! Number of steps in the integral - denominator in the evaluation of the integral
      real*8 funl ! Function subtitution for the first integral
      real*8 funu ! Function subtitution for the second integral
      real*8 sum  ! Intermediate sum in the evaluation of the integrals
      real*8 func ! External function 

      EXTERNAL func

      data pi/3.141592653589793/

!     function subtitution
      funl(x) = 2.d0*x*func(aa+x**2,bb,aa)
      funu(x) = 2.d0*x*func(bb-x**2,bb,aa)
 
!     Zero to variables
      l=0
      outres=0d0
      if(epl.lt.1d-5) goto 333

!     Big loop over 2 segments of the orbit because of the true anomaly
      do ik=1,2

!      Pericenter, apocenter and lower and upper integration limits
       aa=apl*(1d0-epl)
       bb=apl*(1d0+epl)
       b=sqrt(apl*epl)
       a=0.d0

!      Zeros to intermediate sums
       s1=0d0
       s2=0d0

!      Loop over the consencutive iterations of the integral
!      FIRST INTEGRAL
       do n=1,nmax
!       The first evaluation, the same as for the circular orbit of the target
        if (n.eq.1) then
        x=0.5*(a+b)
!       Add 1 to the counter
        l=l+1
        CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,aa+x**2,prob,prob1,theta,phi,ik,outh,outk,outv)      
        outR(l)=aa+x**2

!       Output flush
        do k=1,8 
         outhh(l,k)=outh(k)
         outkk(l,k)=outk(k) 
         outP(l,k)=prob(k)*funl(x)*(b-a)/int(3d0**(nmax-1))/2d0 !Normalization of indidividual impact probabilities
         outTH(l,k)=theta(k) 
         outPH(l,k)=phi(k)
	 outVEL(l,k)=outv(k)
        enddo
        s1=(b-a)*funl(x)*prob1
!      Here comes the integral
       else
         it=int(3d0**(n-2))
        tnm=it
        del=(b-a)/(3*tnm)
        ddel=del+del
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
!       Add 1 to the counter
         l=l+1
         CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,aa+x**2,prob,prob1,theta,phi,ik,outh,outk,outv)
         outR(l)=aa+x**2
!        Output flush
         do k=1,8 
          outhh(l,k)=outh(k)
          outkk(l,k)=outk(k)
          outP(l,k)=prob(k)*funl(x)*(b-a)/2d0/int(3d0**(nmax-1)) !Normalization of indidividual impact probabilities
          outTH(l,k)=theta(k) 
          outPH(l,k)=phi(k)
	  outVEL(l,k)=outv(k)
         enddo
         sum=sum+funl(x)*prob1
         x=x+ddel

!       Add 1 to the counter
         l=l+1
         CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,aa+x**2,prob,prob1,theta,phi,ik,outh,outk,outv)
         outR(l)=aa+x**2

!      Output flush
         do k=1,8 
          outhh(l,k)=outh(k)
          outkk(l,k)=outk(k)
          outP(l,k)=prob(k)*funl(x)*(b-a)/2d0/int(3d0**(nmax-1)) !Normalization of indidividual impact probabilities
          outTH(l,k)=theta(k) 
          outPH(l,k)=phi(k)
	  outVEL(l,k)=outv(k)
         enddo

         sum=sum+funl(x)*prob1
         x=x+del
!        Loop end
11       continue
         s1=(s1+(b-a)*sum/tnm)/3.d0
        endif
       enddo

!      Loop over the consencutive iterations of the integral
!      SECOND INTEGRAL
       do n=1,nmax
!      The first evaluation, the same as for the circular orbit of the target
        if (n.eq.1) then

         x=0.5*(a+b)
!       Add 1 to the counter
         l=l+1
         CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,bb-x**2,prob,prob1,theta,phi,ik,outh,outk,outv)      
         outR(l)=bb-x**2

!        Output flush
         do k=1,8 
          outhh(l,k)=outh(k)
          outkk(l,k)=outk(k)
          outP(l,k)=prob(k)*funu(x)*(b-a)/int(3d0**(nmax-1))/2d0 !Normalization of indidividual impact probabilities
!          outp(l,k)=0d0
         outTH(l,k)=theta(k) 
         outPH(l,k)=phi(k)
	 outVEL(l,k)=outv(k)
        enddo
!         prob1=0d0
        s2=(b-a)*funu(x)*prob1
        else
         it=int(3d0**(n-2))
         tnm=it
         del=(b-a)/(3.*tnm)
         ddel=del+del
         x=a+0.5*del
         sum=0.d0
         do 12 j=1,it
!       Add 1 to the counter
          l=l+1
          CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,bb-x**2,prob,prob1,theta,phi,ik,outh,outk,outv)
          outR(l)=bb-x**2

!       Output flush
          do k=1,8 
           outhh(l,k)=outh(k)
           outkk(l,k)=outk(k)
           outP(l,k)=prob(k)*funu(x)*(b-a)/2d0/int(3d0**(nmax-1)) !Normalization of indidividual impact probabilities
           outTH(l,k)=theta(k) 
           outPH(l,k)=phi(k)
	   outVEL(l,k)=outv(k)
          enddo
          sum=sum+funu(x)*prob1
          x=x+ddel
!       Add 1 to the counter
          l=l+1
          CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,bb-x**2,prob,prob1,theta,phi,ik,outh,outk,outv)
          outR(l)=bb-x**2

!       Output flush
          do k=1,8 
           outhh(l,k)=outh(k)
           outkk(l,k)=outk(k)
           outP(l,k)=prob(k)*funu(x)*(b-a)/2d0/int(3d0**(nmax-1)) !Normalization of indidividual impact probabilities
           outTH(l,k)=theta(k) 
           outPH(l,k)=phi(k)
  	   outVEL(l,k)=outv(k)
          enddo
          sum=sum+funu(x)*prob1
          x=x+del
12       continue
         s2=(s2+(b-a)*sum/tnm)/3d0
        endif

       enddo
      
       outres=(s1+s2)+outres

      enddo

      outres=outres/2d0
      goto 789 ! We're done here, go to 789 to end this

! SIMPLE ADDITIONAL OUTPUT EXAMPLE
!       write(*,*) outres
!      do j=1,coef
!       do k=1,8
!        if (outP(j,k).gt.(0d0)) then
!         write(*,*)outP(j,k),outR(j),outTH(j,k),
!     &   outPH(j,k),outhh(j,k),outkk(j,k)
!        endif
!       enddo
!      enddo
! END OF EXAMPLE
      
      
 333     continue
! This is a light version of the code for cases when the target eccentricity is very small so we don't have to do the whole integral over target's orbit
         outres=0d0
         aa=apl
         x=0d0
         epl=0d0
	 l=1
	 outR(l)=apl
	 CALL KOZAIOPIK(apr,e,inc,omega,apl,epl,aa+x**2,prob,prob1,theta,phi,ik,outh,outk,outv)
          do k=1,8 
           outhh(l,k)=outh(k)
           outkk(l,k)=outk(k)
           outP(l,k)=prob(k)
           outTH(l,k)=theta(k) 
           outPH(l,k)=phi(k)
  	   outVEL(l,k)=outv(k)
          enddo
      
 789 continue
        end



! Function \psi from the article Pokorny et al. 2013 - needed for 
! the total probability

      REAL*8 function func(x,a,b)
      real*8 x,a,b
      func = 1d0/x/sqrt((b-x)*(x-a))*sqrt(a*b)/(asin(1d0)*2)
      end function




! This subroutine computes a collision probability for a target on elliptic 
! non-inclined orbit with a projectile on bound heliocentric orbit
! 
! The subroutine itself computes a probability only for one heliocentric
! distance of the target rpl, the total probability must be computed as 
! an integral of individual probabilities from pericenter to apocenter



      SUBROUTINE KOZAIOPIK(a,e,inc,om,apl,epl,rpl,prob,prob1,theta,phi,ik,outh,outk,outv)


      IMPLICIT NONE
!     INPUT VARIABLES
      real*8 a ! Semimajor axis of the projectile [AU]
      real*8 e ! Eccentricity of the projectile 
      real*8 inc ! Inclination of the projectile [rad]
      real*8 om ! Argument of pericenter of the projectile [rad]
      real*8 apl ! Semimajor axis of the target [AU]
      real*8 epl ! Eccentricity of the target 
      real*8 rpl ! Heliocentric distance of the target [AU]
!     END OF INPUT

!     OUTPUT VARIABLES
      real*8 prob(8) ! Intrinsic collisional probabilities for all impact configurations [1D array]
      real*8 prob1 ! Total intrinsic collisional probability = sum of prob(8)
      real*8 outk(8) ! Array of all obtained k=e*cos(omega) in the impact configurations
      real*8 outh(8) ! Array of all obtained h=e*sin(omega) in the impact configurations
      real*8 theta(8) ! Ecliptic latitude of the impact radiants for the impact configurations
      real*8 phi(8) ! Ecliptic longitude of the impact radiants for the impact configurations
      real*8 outv(8) ! Impact velocity for the impact configurations
!     END OF OUTPUT

      complex*16 x(3) ! Roots of the system of equations (10) and (16) = impact configurations - cubic equation

      real*8 Ca ! Coefficients of the cubic Equation
      real*8 Cb
      real*8 Cc
      real*8 Cd

      real*8 y(4) ! Real part of the roots - maximum 8 impact configurations => 4 unique values of k=e*cos(omega) 
      real*8 z(8) ! h=e*cos(omega) - from EQ (10) in Vokrouhlicky, Pokorny, Nesvorny 2012
      real*8 k ! k=e*cos(omega) - just more convenient and readable variable without array
      real*8 hh(2) ! Array of 2 possible values of h=e*sin(omega) for one unique k=e*cos(omega)
      real*8 h ! Selected h from array hh(2)

      real*8 kk ! Initial value of k=e*cos(omega) from the initial conditions
      real*8 hhh ! Initial value of h=e*sin(omega) from the initial conditions

      real*8 CapC ! Integral of the motion EQ (16) in Pokorny, Vokrouhlicky 2013
      real*8 c ! Integral of the motion c = eta * cos (inc)
      real*8 eta ! \eta = \sqrt{1-e^2}
      real*8 gamm ! Factor Gamma for scaling the lenght of the Lidov-Kozai cycle - see below EQ (28) in Vokrouhlicky, Pokorny, Nesvorny 2012

      real*8 dk ! Displacement in k, for an exact impact configuration point - EQ (19) in Vokrouhlicky, Pokorny, Nesvorny 2012
      real*8 dh ! Displacement in h, for an exact impact configuration point - EQ (19) in Vokrouhlicky, Pokorny, Nesvorny 2012
      real*8 dkdt ! Derivation of k=e*cos(omega) by time - EQ (27) in Vokrouhlicky, Pokorny, Nesvorny 2012 
      real*8 dhdt ! Derivation of h=e*sin(omega) by time - EQ (28) in Vokrouhlicky, Pokorny, Nesvorny 2012 

      real*8 phk ! EQ (21) in Vokrouhlicky, Pokorny, Nesvorny 2012 
      real*8 phh ! EQ (22) in Vokrouhlicky, Pokorny, Nesvorny 2012 
      real*8 psk ! EQ (23) in Vokrouhlicky, Pokorny, Nesvorny 2012 
      real*8 psh ! EQ (24) in Vokrouhlicky, Pokorny, Nesvorny 2012 
      real*8 den ! EQ (20) in Vokrouhlicky, Pokorny, Nesvorny 2012 

      real*8 Pom ! T_Kozai - The lenght of the one Lidov-Kozai cycle [yr] - see Kinoshita & Nakai 2007

      real*8 P2 ! EQ (17) in Pokorny, Vokrouhlicky 2013 - Probability part 2

      real*8 aph
      real*8 rho2
      real*8 pi

      real*8 kx
      real*8 hx
      real*8 dkdt1
      real*8 dhdt1
      real*8 dist
      real*8 timm
      real*8 r1
      real*8 r2
      real*8 aa
      real*8 th
      real*8 ph
      real*8 Ie
      real*8 Iinc
      real*8 Iom
      real*8 vel

      integer i
      integer j
      integer l
      integer m
      integer m2
      integer nroot
      integer znk(2)
      integer ddd
      integer CUBIC
      integer ik
      data pi/3.141592653589793/

      real*8 err
      real*8 preci

!     Pericenter r1 and apocenter r2 of the target
      r1=apl*(1d0-epl)
      r2=apl*(1d0+epl)

!     Zeros to variables
      theta=0d0
      phi=0d0      
      prob=0d0
      prob1=0d0
      outv=0d0
      ddd=0

!     Determine initial k,h for the projectile
      kk=e*cos(om)
      hhh=e*sin(om)
      
!     Determine aplha,eta
      aph=rpl/a
      eta=sqrt(1d0-e*e)

!     Evaluate integrals of motion       
      c=eta*cos(inc)
      CapC=1./eta/eta*((2d0+3.*e*e)*(3d0*c*c-eta*eta)+15d0*(eta*eta-c*c)*(kk*kk-hhh*hhh))

!     Gamma factor, where 0.0009546 is mass of the Jupiter in Solar Masses
!     and 140.954423151 - this value is not important for the collision probability.
!     However scales the length of Lidov-Kozai cycle
      gamm=2d0*pi*0.0009546*sqrt(a**3)/140.954423151     

!     Iteration value      
      l=1

!!!   Should be computed only once (CHANGE)
      CALL KINOSHITA(CapC,c,eta,om,gamm,pom)

!     Ascending/Descending node i = 1 Ascending, i = 2 Descending
 
      do i=1,2

!      Determine coefficients for the Cubic equation

       Ca=(-1.)**(i+1)*30.*aph
       Cb=6.*(3.*aph*aph+5.*aph-5.*c*c)
       Cc=(-1.)**(i+1)*aph*(36.*aph-20.-24.*c*c-CapC)
       Cd=18.*aph*aph-20*aph-24*aph*c*c+30*c*c-CapC*aph      

!      Call a routine to solve the cubic equation
       nroot = CUBIC(Ca, Cb, Cc, Cd, x)

!       Loop over the roots of the cubic equation ie. the impact configurations
        do j=1, nroot

!        We want only non-complex solutions 
         if ((aimag(x(j)).eq.(0.)).and.(abs(real(x(j))).lt.(1.))) then 
!        Get the real part of the complex root of the equation, k=e*cos(omega)
         y(l)= real(x(j))
!        Evaluate h=e*cos(omega) - EQ (10) in Vokrouhlicky, Pokorny, Nesvorny 2012
         z(2*l-1)=sqrt(1-y(l)**2-aph*(1+(-1.)**(i+1)*y(l)))
         z(2*l-0)=-z(2*l-1)

!        Check if h=e*cos(omega) can exist
         if ((1-y(l)**2-aph*(1+(-1.)**(i+1)*y(l))).lt.(0d0)) then
          goto 200 ! GOTO END OF THE MAIN LOOP (it should never happen, however :))
         endif


!        Assign to more conveniet(readable) variables
         k=y(l)
         hh(1)=z(2*l-1)
         hh(2)=z(2*l-0)
!        Set number of unique h=e*cos(omega)
         m2=2 
!        Look whether h=e*cos(omega) is unique, if not m2=1
         if (hh(2).eq.hh(1)) m2=1



!        Forced oscillations - max 4 impact configurations
         if (CapC.lt.(2*(3*c**2-1))) then
          m2=1
          hh(1)=abs(hh(1))
!         See, whether omega is >180°, if yes then h must be < 0 
          if (om.gt.pi) hh(1)=-hh(1)
         endif

!        Loop over the impact configurations
         do m=1,m2
!        Even more readable variables
         h=hh(m)
         eta=sqrt(1-h**2-k**2)

!        Impact eccentricity,inclination and argument of pericenter
         Ie=sqrt(k**2+h**2)
         Iinc =acos(c/eta)
         Iom = acos(k/Ie)
!        Factor /beta from Eq (21) - see below the equation for a formula
         aa=(-1.)**(ik+1)*sqrt((rpl-r1)*(r2-rpl))/(apl*sqrt(1d0-epl**2))

!        Solution for the ambiquity of acos in the determination of argument of pericenter of the impact
         if (h.lt.(0.0)) Iom = 2*pi - Iom

!        Call Radiant and P2 determination (EQ. 17) - computes also impact velocities, which may be also interesting output
         CALL RADIANT_ECC(a,Ie,Iinc,Iom,apl,epl,rpl,th,ph,P2,ik,vel)
!        Outputs from radiants and also k,h for intersections to arrays
         theta(l*2-2+m)=th
         phi(l*2-2+m)=ph
         outh(l*2-2+m)=h
         outk(l*2-2+m)=k
         outv(l*2-2+m)=vel
!        Factor /Lambda from EQ 23 (formerly named /rho_2 thus may be misleading (TO BE CHANGED FUTURE)
         rho2=sqrt( &
     &       (eta**2 - c**2) / &
     &       ( &
     &       (1d0+aa**2)*(eta**2-c**2)*(1+(-1.)**(i+1)*k)**2 &
     &     + ( &
     &       (h*eta)+(-1.)**(i+1)*aa*c*(1+(-1.)**(i+1)*k) &
     &       )**2 &
     &       ) &
     &       )
!         EQ 22
          dist=(eta**2-aph*(1+(-1.)**(i+1)*k))*rho2
!         Partial derivations of the distance EQs (21) - (24) in Vokrouhlicky, Pokorny, Nesvorny 2012 (ICARUS)
          phk=2.*k*(2.*(7.-3.*c*c+3.*h*h-12.*k*k)+CapC)
          phh=2.*h*(2.*(-8.+12.*c*c+3.*k*k+18.*h*h)+CapC)
          psk=-2.*(k+(-1.)**(i+1)*aph/2.)*rho2
          psh=-2.*(h)*rho2
!         Denominator EQ (20) in Vokrouhlicky, Pokorny, Nesvorny 2012 (ICARUS) 
          den=(phh*psk-phk*psh)
!         And finally displacements EQ (19) in Vokrouhlicky, Pokorny, Nesvorny 2012 (ICARUS)
          dk=phh/den/a
          dh=-phk/den/a

          kx=k
          hx=h

          timm=0

!     If Denominator is too small we have to compute the time in the proximity of planet numerically (it is a bit slower, but we have quite large steps)
!       if (abs(den).lt.(2E-3)) then
      if (abs(den).lt.(2E-3)) then
!      We do not know the radius of the target, thus we set this value as an approximate value (but should be ok for most of cases).      
       err=1D-9
       preci=1D-6
!      Set zero values for a special case commented in Vokrouhlicky, Pokorny, Nesvorny 2012 (ICARUS) - Appendix
       znk(1)=0
       znk(2)=0


       do while (abs(dist).lt.(err))
        dkdt1=-3d0/2d0*preci*gamm*eta*h*(1- &
     &  5d0/2d0*k*k*(eta*eta-c*c)/(1-eta*eta)/eta/eta+  &
     &  5d0/2d0*(c*c-eta**4d0)/(eta**4d0)*h*h/(1d0-eta*eta))
        dhdt1=3d0/2d0*preci*gamm*eta*k*(1+5d0/2d0*h*h*c*c/eta**4d0)      
        k=k-dkdt1
        h=h-dhdt1
        eta=sqrt(1-h*h-k*k)
        dist=(eta*eta-aph*(1+(-1.)**(i+1)*k))* &
     &        sqrt((eta**2 - c**2) / &
     &       ( &
     &       (1d0+aa**2)*(eta**2-c**2)*(1+(-1.)**(i+1)*k)**2 &
     &     + (  &
     &       (h*eta)+(-1.)**(i+1)*aa*c*(1+(-1.)**(i+1)*k)  &
     &       )**2  &
     &       ))
!       Dist should be 
        if (dist.ge.(0.0)) znk(1)=1
        if (dist.le.(0.0)) znk(2)=1
        timm=timm+preci
             enddo
        if ((znk(1).eq.(1)).and.(znk(2).eq.(1))) ddd=1
        znk(1)=0
        znk(2)=0
!        write(*,*) timm, "END 1", ddd
        k=kx
        h=hx
        eta=sqrt(1-h*h-k*k)
        dist=(eta*eta-aph*(1+(-1.)**(i+1)*k))*  &
     &    sqrt((eta**2 - c**2) /  &
     &  (  &
     &  (1d0+aa**2)*(eta**2-c**2)*(1+(-1.)**(i+1)*k)**2  &
     &  + (  &
     &     (h*eta)+(-1.)**(i+1)*aa*c*(1+(-1.)**(i+1)*k)  &
     &     )**2  &
     &  ))
!        write(*,*) dist
            do while (abs(dist).lt.(err))
        dkdt1=-3./2.*preci*gamm*eta*h*(1-  &
     &    5./2.*k*k*(eta*eta-c*c)/(1-eta*eta)/eta/eta+  &
     &    5./2.*(c*c-eta**4.)/(eta**4.)*h*h/(1-eta*eta))
        dhdt1=3./2.*preci*gamm*eta*k*(1+5./2.*h*h*c*c/eta**4.)      
        k=k+dkdt1
        h=h+dhdt1
        eta=sqrt(1-h*h-k*k)
        dist=(eta*eta-aph*(1+(-1.)**(i+1)*k))*  &
     &    sqrt((eta**2 - c**2) /  &
     &  (  &
     &  (1d0+aa**2)*(eta**2-c**2)*(1+(-1.)**(i+1)*k)**2  &
     &  + (  &
     &     (h*eta)+(-1.)**(i+1)*aa*c*(1+(-1.)**(i+1)*k)  &
     &     )**2  &
     &  ))
        if (dist.ge.(0.0)) znk(1)=1
        if (dist.le.(0.0)) znk(2)=1
        timm=timm+preci


! !         write(*,*) dist,k,h,timm
! 	call sleep(0)
             enddo

        if ((znk(1).eq.(1)).and.(znk(2).eq.(1))) ddd=1
        znk(1)=0
        znk(2)=0
!        write(*,*) timm, "END 2",ddd
        k=kx
        h=hx
        if (ddd.eq.(1)) then
        timm=timm/2d0
!       write(*,*) "DOUBLE MODE", timm,abs(timm)/pom,P2
        else
!        write(*,*) "SINGLE MODE", timm,abs(timm)/pom,P2
        endif
        timm=timm/err

        prob(l*2-2+m) =    &
     &    abs(timm)/pom*P2/a

!!!! 1/a is from R/a = rho

!     Denominator is OK, so we can go straighforward as in the article
      else
        k=kx
        h=hx

        dkdt=-3./2.*gamm*eta*h*(1-  &
     &    5./2.*k*k*(eta*eta-c*c)/(1-eta*eta)/eta/eta+  &
     &    5./2.*(c*c-eta**4.)/(eta**4.)*h*h/(1-eta*eta))
        dhdt=3./2.*gamm*eta*k*(1+5./2.*h*h*c*c/eta**4.)  
            prob(l*2-2+m)=abs(dk/dkdt)*2/pom*P2
      endif
      prob1=prob(l*2-2+m)+prob1
      enddo


              l=l+1
          endif
 200  continue
        enddo
      

      enddo

      END

     
      SUBROUTINE RADIANT_ECC(a,e,inc,om,apl,epl,rpl,theta,phi,P2,k,vel)
      IMPLICIT NONE

! This subroutine computes radiants of a particle (phi,theta)
! for a eccentric coplanar target and projectile
! apex of the target's motion has theta = 0, phi = 0
! theta means ecliptic latitude 
! phi means ecliptic longitude      
!!! INPUT/OUTPUT



! INPUT Variables: 
! a ... semimajor axis of the projectile
! e ... eccentricity of the projectile
! inc ... inclination of the projectile in radians
! om ... argument of pericenter of the projectile in radians
! apl ... semimajor axis of the target
! epl ... eccentricity of the target
! rpl ... heliocentric distance of the target

      real*8 a,e,inc,om,apl,epl,rpl
      integer k

! OUTPUT Variables:
! theta ... ecliptic latitude (-pi/2,pi/2)
! phi ..... ecliptic longitude (-pi,pi)
! P2 ...... second term in probability calculation

      real*8 theta,phi,P2


! PROCESSING VARIABLES

      real*8 eta,etapl,asc,dsc,GM,r1,r2
      real*8 v1,v2,v3,vel,pi,rad,bigR,bigF,P,R2F2

!!! DOUBLE PRECISION
      data pi/3.141592653589793/

! GM = gravitational constat * mass of the gravitational center
! in [astronomical unit**3 * mean solar day ** (-2) * solar mass**(-1) 
      GM = 0.2959122104742908d-03       
! Some initial computations
      rad=pi/180.0d0
      eta=sqrt(1d0-e*e)
      etapl=sqrt(1d0-epl**2)
      r1=apl*(1d0-epl)
      r2=apl*(1d0+epl)
      bigR=(-1.)**(k+1)*sqrt((rpl-r1)*(r2-rpl))/rpl
      bigF=etapl*apl/rpl
      P=a*eta**2/rpl
      R2F2=bigR**2+bigF**2

      asc = P/(1d0+e*cos(om))
      dsc = P/(1d0-e*cos(om))

!     Zeros to variables
      v1=0d0
      v2=0d0
      v3=0d0
      vel=0d0

      if ((1d0-1e-4.le.asc).and.(asc.le.1d0+1e-4)) then

      v1 = -(sqrt(apl/rpl*P)*(bigF*cos(inc) - e*sin(om)/P*bigR)  &
     &      - R2F2) / sqrt(R2F2)
      v2 = (e*sin(om)*bigF/P + bigR*cos(inc)) * sqrt(apl/rpl*P)  &
     &      / sqrt(R2F2)
      v3 = - sqrt(apl/rpl*P) * sin(inc)

      vel=sqrt(v1**2+v2**2+v3**2)
      theta=asin(v3/vel)
      phi=atan2(v2,v1)
!      write(*,*)theta/rad,phi/rad, "asc"
      endif


      if ((1d0-1e-4.le.dsc).and.(dsc.le.1d0+1e-4)) then

      v1 = -(sqrt(apl/rpl*P)*(bigF*cos(inc) + e*sin(om)/P*bigR)  &
     &      - R2F2) / sqrt(R2F2)
      v2 = ( - e*sin(om)*bigF/P + bigR*cos(inc)) * sqrt(apl/rpl*P)  &
     &      / sqrt(R2F2)
      v3 = + sqrt(apl/rpl*P) * sin(inc)


      vel=sqrt(v1**2+v2**2+v3**2)
      theta=asin(v3/vel)
      phi=atan2(v2,v1)
      endif
      
      theta=theta/rad
      phi=phi/rad
      P2=1d0/apl/4d0/sqrt(R2F2)*vel/sqrt(vel**2-v1**2)/sqrt(a**3)
      vel=vel*(apl*dsqrt(GM)*149597870.7/86400.0)
      END

!       ==================================================

        SUBROUTINE COMELP(HK,CK,CE)

!       ==================================================
!       Purpose: Compute complete elliptic integrals K(k)
!                and E(k)
!       Input  : K  --- Modulus k ( 0 ?? k ?? 1 )
!       Output : CK --- K(k)
!                CE --- E(k)
!       ==================================================
!
        IMPLICIT NONE
        real*8 PK,HK,CK,CE,AE,BE,AK,BK
        
        PK=1.0-HK*HK
        IF (HK.EQ.1.0) THEN
           CK=1.0+300.
           CE=1.0
        ELSE
           AK=(((.01451196212*PK+.03742563713)*PK  &
     &        +.03590092383)*PK+.09666344259)*PK+  &
     &        1.38629436112
           BK=(((.00441787012*PK+.03328355346)*PK+  &
     &        .06880248576)*PK+.12498593597)*PK+.5
           CK=AK-BK*dLOG(PK)
           AE=(((.01736506451*PK+.04757383546)*PK+  &
     &        .0626060122)*PK+.44325141463)*PK+1.0
           BE=(((.00526449639*PK+.04069697526)*PK+  &
     &        .09200180037)*PK+.2499836831)*PK
           CE=AE-BE*LOG(PK)
           CK=AK-BK*dlog(PK)
        ENDIF

        RETURN
        END

!       ==================================================

      SUBROUTINE KINOSHITA(CapC,c,eta,om,gamm,pom)


! Kinoshita-Nakai 2007
!
! Input variables
! CapC .. energy integral in quadrupole Kozai oscillations
! c ..... particle's integral of motion c = eta * cos(i)
! eta ... eta = (1 - e**2)
! om .... argument of pericenter of the particle
!
!          Ms/Mj = 0.0009546, aJ**3 = 140.954423151 AU**3
! gamm ... gamm=2d0*pi*0.0009546*sqrt(a*a*a)/140.954423151
!
      implicit none
      real*8 CapC,c,eta,om,gamm
! Output variables
! pom .... one period of argument of pericenter in years
!          (if in forced libration mode we take only 1/2 of it)

      real*8 pom

! Temp variables

      real*8 C1,C2,x0,x1,x2,DD,tmp,a0,a1,a2,HK,CK,CE,nom,pi
      complex*16 xq(2)
        data pi/3.141592653589793/

        a0=0d0
        a1=0d0
        a2=0d0
        x1=0d0
        x2=0d0
        
      C1=5d0+5d0*c*c
      C2=5d0*c*c/eta/eta+eta*eta+5d0*(1-eta*eta)*(1-c*c/eta/eta)*  &
     &  cos(2*om)
      x0=1d0/4d0*(C1-C2)

      DD = 1d0/2d0*(C1+C2)*1d0/2d0*(C1+C2)-12d0*5d0*c*c
      if(DD .ge. 0.)then
          xq(1)=cmplx((-1d0/2d0*(C1+C2)+sqrt(DD))/2./(-3d0), 0.)
          xq(2)=cmplx((-1d0/2d0*(C1+C2)-sqrt(DD))/2./(-3d0), 0.)
        else
          xq(1)=cmplx(-1d0/2d0*(C1+C2)/2./(-3d0), +sqrt(-DD)/2./(-3d0))
          xq(2)=cmplx(-1d0/2d0*(C1+C2)/2./(-3d0), -sqrt(-DD)/2./(-3d0))
        endif
            if (aimag(xq(1)).eq.(0.)) x1=real(xq(1))
            if (aimag(xq(2)).eq.(0.)) x2=real(xq(2))
      if (x1.gt.x2) then
      tmp=x1
      x1=x2
      x2=tmp
      endif
      if (x0.le.x1) then
      a0=x0
      a1=x1
      a2=x2

      endif
      if (x0.ge.x2) then
      a0=x1
      a1=x2
      a2=x0
      endif
      if ((x0.ge.x1).and.(x0.le.x2)) then
      a0=x1
      a1=x0
      a2=x2
      endif

!   Complete Elliptic Integral of the First Kind = CK
      HK=sqrt((a1-a0)/(a2-a0))

        CALL COMELP(HK,CK,CE)
!      write(*,*) CK
!       IF (HK.NE.1.0) WRITE(*,*) HK,CK,CE
!       IF (HK.EQ.1.0) WRITE(*,*) HK,CE
      Nom=3d0*sqrt(6d0)*pi/8d0/CK*sqrt(a2-a0)*gamm
      Pom=2d0*pi/nom


!      Variable Pom is a period of one Kozai cycle (see Kinoshita-Nakai 2007)

      if (CapC.lt.(2*(3*c*c-1))) pom=pom/2d0


      END



      FUNCTION quad(a, b, c, x)

! Solve a quadratic equation where a, b, and c are real.
!   a*x*x + b*x + c = 0
! This function returns the number of roots
!
! Public Variables
!   a, b, c     ... coefficients (input)
!   x(i)        ... two complex*16       solutions (output)
!   nroot       ... number of roots (output)
! ----------------------------------------------------------------------

       IMPLICIT NONE
      complex*16 x(2)
      real*8 a,b,c,DD
      integer QUAD,NROOT
      
      if(a .eq. 0.)then
        if(b .eq. 0.)then
!           We have a non-equation; therefore, we have no valid solution
          nroot = 0
        else
!           We have a linear equation with 1 root.
          nroot = 1
          x(1) = cmplx(-c/b, 0.)
        endif
      else
!       We have a true quadratic equation.  Apply the quadratic formula to find two roots.
      nroot = 2
        DD = b*b-4.*a*c
      if(DD .ge. 0.)then
          x(1) = cmplx((-b+sqrt(DD))/2./a, 0.)
          x(2) = cmplx((-b-sqrt(DD))/2./a, 0.)
        else
          x(1) = cmplx(-b/2./a, +sqrt(-DD)/2./a)
          x(2) = cmplx(-b/2./a, -sqrt(-DD)/2./a)
        endif
      endif

! Return the number of roots
      QUAD = NROOT

      end

! ----------------------------------------------------------------------
      FUNCTION cubic(a, b, c, d, x)
! ----------------------------------------------------------------------
!   Solve a cubic equation where a, b, c, and d are real.
!     a*x**3 + b*x**2 + c*x + d = 0
!   This function returns the number of roots
!  
!   Public Variables
!     a, b, c, d  ... coefficients (input)
!     x(i)        ... three (generally) complex*16 solutions (output)
!     nroot       ... number of roots (output)
!   Local Variables:
!     y1, y2, y3  ... three transformed solutions
!  
!   Formula used are given in Tuma, "Engineering Mathematics Handbook", p7
!     (McGraw Hill, 1978).
!     Step 0: If a is 0. use the quadratic formula to avoid dividing by 0.
!     Step 1: Calculate p and q
!             p = ( 3*c/a - (b/a)**2 ) / 3
!             q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
!     Step 2: Calculate discriminant D
!             D = (p/3)**3 + (q/2)**2
!     Step 3: Depending on the sign of D, we follow different strategy.
!             If D<0, three distinct real roots.
!             If D=0, three real roots of which at least two are equal.
!             If D>0, one real and two complex*16 roots.
!     Step 3a: For D>0 and D=0,
!             Calculate u and v
!             u = cubic_root(-q/2 + sqrt(D))
!             v = cubic_root(-q/2 - sqrt(D))
!             Find the three transformed roots
! 
!             y1 = u + v
!             y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
!             y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
!     Step 3b Alternately, for D<0, a trigonometric formulation is more convenient
!             y1 =  2 * sqrt(|p|/3) * cos(phi/3)
!             y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
!             y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
!             where phi = acos(-q/2/sqrt(|p|**3/27))
!                   pi  = 3.141592654...
!     Step 4  Finally, find the three roots
!             x = y - b/a/3
!   ----------------------------------------------------------------------

       IMPLICIT NONE

!   Declare variables
      complex*16 x(3)

      real*8 DD,q,p,temp1,phi,temp2,u,v,y1,y2,y3,y2r,y2i,a,b,c,d,pi

      integer QUAD,NROOT,CUBIC


      data pi/3.141592653589793/
!   Step 0: If a is 0 use the quadratic formula. -------------------------
      if(a .eq. 0.)then
        nroot = QUAD(b, c, d, x)
        CUBIC = NROOT
        return
      endif

!   Cubic equation with 3 roots
      nroot = 3

!   Step 1: Calculate p and q --------------------------------------------
      p  = c/a - b*b/a/a/3.
      q  = (2.*b*b*b/a/a/a - 9.*b*c/a/a + 27.*d/a) / 27.

!   Step 2: Calculate DD (discriminant) ----------------------------------
      DD = p*p*p/27. + q*q/4.

!   Step 3: Branch to different algorithms based on DD -------------------
      if(DD .lt. 0.)then
!         Step 3b:
!         3 real unequal roots -- use the trigonometric formulation
        phi = acos(-q/2./sqrt(abs(p*p*p)/27.))
        temp1=2.*sqrt(abs(p)/3.)
        y1 =  temp1*cos(phi/3.)
        y2 = -temp1*cos((phi+pi)/3.)
        y3 = -temp1*cos((phi-pi)/3.)
      else
!         Step 3a:
!         1 real root & 2 conjugate complex*16 roots OR 3 real roots (some are equal)
        temp1 = -q/2. + sqrt(DD)
        temp2 = -q/2. - sqrt(DD)
        u = exp(log(abs(temp1))/3.)
        v = exp(log(abs(temp2))/3.)
        if(temp1 .lt. 0.) u=-u
        if(temp2 .lt. 0.) v=-v
        y1  = u + v
        y2r = -(u+v)/2.
        y2i =  (u-v)*sqrt(3.)/2.
      endif

!   Step 4: Final transformation -----------------------------------------
      temp1 = b/a/3.
      y1 = y1-temp1
      y1=y1-(d+y1*(c+y1*(b+y1*a)))  &
     &  /(c+y1*(2.*b+y1*3.*a))
      y2 = y2-temp1
      y2=y2-(d+y2*(c+y2*(b+y2*a)))  &
     &  /(c+y2*(2.*b+y2*3.*a))
      y3 = y3-temp1
      y3=y3-(d+y3*(c+y3*(b+y3*a)))  &
     &  /(c+y3*(2.*b+y3*3.*a))
      y2r=y2r-temp1


!   Assign answers -------------------------------------------------------
      if(DD .lt. 0.)then
        x(1) = cmplx( y1,  0.)
        x(2) = cmplx( y2,  0.)
        x(3) = cmplx( y3,  0.)
      elseif(DD .eq. 0.)then
        x(1) = cmplx( y1,  0.)
      y2r=y2r-(d+y2r*(c+y2r*(b+y2r*a)))  &
     &  /(c+y2r*(2.*b+y2r*3.*a))
        x(2) = cmplx(y2r,  0.)
        x(3) = cmplx(y2r,  0.)
      else
        x(1) = cmplx( y1,  0.)
        x(2) = cmplx(y2r, y2i)
        x(3) = cmplx(y2r,-y2i)
      endif

!   Return the number of roots
      CUBIC = NROOT

      end


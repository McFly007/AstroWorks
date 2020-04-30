  module my_functions
 contains

 FUNCTION PHI1(X)
  REAL*8  :: PHI1
  REAL*8, INTENT(IN) :: X
  INTEGER :: PICK
  REAL*8, DIMENSION(7) :: X1,X2,Y1,Y2,D1,D2
  REAL*8 T,A,B,DEG2RAD
  DEG2RAD=0.01745329251d0

  if(X.le.7.5d0)  PICK=1
  if((x.gt.7.5).and.(x.le.30.0))  PICK=2
  if((x.gt.30.0).and.(x.le.60.0))  PICK=3
  if((x.gt.60.0).and.(x.le.90.0))  PICK=4
  if((x.gt.90.0).and.(x.le.120.0))  PICK=5
  if((x.gt.120.0).and.(x.le.150.0))  PICK=6
  if((x.gt.150.0).and.(x.le.180.0))  PICK=7
  
  X1=(/0d0,7.5d0,30d0,60d0,90d0,120d0,150d0/)
  X2=(/7.5d0,30d0,60d0,90d0,120d0,150d0,180d0/)
  Y1=(/1d0,7.5d-1,3.3486016d-1,1.3410560d-1,5.1104756d-2,2.1465687d-2,3.6396989d-3/)
  Y2=(/7.5d-1,3.3486016d-1,1.3410560d-1,5.1104756d-2,2.1465687d-2,3.6396989d-3,0d0/)
  D1=(/-1.9098593d0,-1.9098593d0,-5.5463432d-1,-2.4404599d-1,-9.4980438d-2,-2.1411424d-2,-9.1328612d-2/)
  D2=(/-1.9098593d0,-5.5463432d-1,-2.4404599d-1,-9.4980438d-2,-2.1411424d-2,-9.1328612d-2,0d0/)
  
  T=(X*DEG2RAD-X1(PICK)*DEG2RAD)/(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)
  A=D1(PICK)*(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)-(Y2(PICK)-Y1(PICK))
  B=-D2(PICK)*(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)+(Y2(PICK)-Y1(PICK))
  PHI1=(1d0-T)*Y1(PICK)+T*Y2(PICK)+T*(1d0-T)*((1d0-T)*A+B*T)
 END FUNCTION PHI1

  FUNCTION PHI2(X)
  REAL*8  :: PHI2
  REAL*8, INTENT(IN) :: X
  INTEGER :: PICK
  REAL*8, DIMENSION(7) :: X1,X2,Y1,Y2,D1,D2
  REAL*8 T,A,B,DEG2RAD
  DEG2RAD=0.01745329251d0

  if(X.le.7.5d0)  PICK=1
  if((x.gt.7.5).and.(x.le.30.0))  PICK=2
  if((x.gt.30.0).and.(x.le.60.0))  PICK=3
  if((x.gt.60.0).and.(x.le.90.0))  PICK=4
  if((x.gt.90.0).and.(x.le.120.0))  PICK=5
  if((x.gt.120.0).and.(x.le.150.0))  PICK=6
  if((x.gt.150.0).and.(x.le.180.0))  PICK=7
  
  X1=(/0d0,7.5d0,30d0,60d0,90d0,120d0,150d0/)
  X2=(/7.5d0,30d0,60d0,90d0,120d0,150d0,180d0/)
  Y1=(/1d0,9.25d-1,6.2884169d-1,3.1755495d-1,1.2716367d-1,2.2373903d-2,1.6505689d-4/)
  Y2=(/9.25d-1,6.2884169d-1,3.1755495d-1,1.2716367d-1,2.2373903d-2,1.6505689d-4,0d0/)
  D1=(/-5.7295780d-1,-5.7295780d-1,-7.6705367d-1,-4.5665789d-1,-2.8071809d-1,-1.1173257d-1,-8.6573138d-8/)
  D2=(/-5.7295780d-1,-7.6705367d-1,-4.5665789d-1,-2.8071809d-1,-1.1173257d-1,-8.6573138d-8,0d0/)
  
  T=(X*DEG2RAD-X1(PICK)*DEG2RAD)/(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)
  A=D1(PICK)*(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)-(Y2(PICK)-Y1(PICK))
  B=-D2(PICK)*(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)+(Y2(PICK)-Y1(PICK))
  PHI2=(1d0-T)*Y1(PICK)+T*Y2(PICK)+T*(1d0-T)*((1d0-T)*A+B*T)
 END FUNCTION PHI2
 
   FUNCTION PHI3(X)
  REAL*8  :: PHI3
  REAL*8, INTENT(IN) :: X
  INTEGER :: PICK
  REAL*8, DIMENSION(8) :: X1,X2,Y1,Y2,D1,D2
  REAL*8 T,A,B,DEG2RAD
  DEG2RAD=0.01745329251d0

  if(X.le.0.3d0)  PICK=1
  if((x.gt.0.3).and.(x.le.1d0))  PICK=2
  if((x.gt.1.0).and.(x.le.2d0))  PICK=3
  if((x.gt.2.0).and.(x.le.4d0))  PICK=4
  if((x.gt.4.0).and.(x.le.8d0))  PICK=5
  if((x.gt.8.0).and.(x.le.12d0))  PICK=6
  if((x.gt.12.0).and.(x.le.20d0))  PICK=7
  if((x.gt.20.0).and.(x.le.180d0))  PICK=8
  
  X1=(/0d0,0.3d0,1d0,2d0,4d0,8d0,12d0,20d0/)
  X2=(/0.3d0,1d0,2d0,4d0,8d0,12d0,20d0,30d0/)
  Y1=(/1d0,8.3381185d-1,5.7735424d-1,4.2144772d-1,2.3174230d-1,1.0348178d-1,6.1733473d-2,1.6107006d-2/)
  Y2=(/8.3381185d-1,5.7735424d-1,4.2144772d-1,2.3174230d-1,1.0348178d-1,6.1733473d-2,1.6107006d-2,0d0/)
  D1=(/-1.0630097d-1,-4.1180439d1,-1.0366915d1,-7.5784615d0,-3.6960950d0,-7.8605652d-1,-4.6527012d-1,-2.0459545d-1/)
  D2=(/-4.1180439d1,-1.0366915d1,-7.5784615d0,-3.6960950d0,-7.8605652d-1,-4.6527012d-1,-2.0459545d-1,0d0/)
  

  T=(X*DEG2RAD-X1(PICK)*DEG2RAD)/(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)
  A=D1(PICK)*(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)-(Y2(PICK)-Y1(PICK))
  B=-D2(PICK)*(X2(PICK)*DEG2RAD-X1(PICK)*DEG2RAD)+(Y2(PICK)-Y1(PICK))
  PHI3=(1d0-T)*Y1(PICK)+T*Y2(PICK)+T*(1d0-T)*((1d0-T)*A+B*T)
  if(X.gt.30) PHI3=0d0
 END FUNCTION PHI3
 
 FUNCTION scalar(a, b)
  REAL*8 :: scalar
  REAL*8, DIMENSION(3), INTENT(IN) :: a, b

  scalar = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 END FUNCTION scalar

 end module my_functions
 
 PROGRAM ASTEROID_TYPES
 use my_functions
 implicit none
 
 real*8 t,r,delta,alpha,lambda,beta
 real*8 pi,VIS1AU,VISUAL
 real*8 DIA(5),G1(6),G2(6),GEOA,H(5,6),VISLINE(30)
 integer i,j,k1,k2,k3
 
 
 pi=dacos(-1d0)
 G1(1)=0.26
 G1(2)=0.27
 G1(3)=0.15
 G1(4)=0.82
 G1(5)=0.83
 G1(6)=0.96
 
 G2(1)=0.38
 G2(2)=0.35
 G2(3)=0.60
 G2(4)=0.02
 G2(5)=0.05
 G2(6)=0.02
 
 
 ! Values from Shevchenko et al. (2016)
 do k1=1,6
 if(k1.eq.1) GEOA=0.22
 if(k1.eq.2) GEOA=0.17
 if(k1.eq.3) GEOA=0.45
 if(k1.eq.4) GEOA=0.061
 if(k1.eq.5) GEOA=0.042
 if(k1.eq.6) GEOA=0.049
 do k2=1,5
 DIA(k2)=k2*0.5d0
 H(k2,k1)=-5d0*log10(DIA(K2)*dsqrt(GEOA)/1329d0)
 enddo
 enddo
 
 open(243,file="Datafile.dat", status="OLD")
 open(244,file="Visual_Magnitudes", status="UNKNOWN")
 
 read(243,*)
 write(244,*) "#VISUAL MAGNITUDES 5 diameters for each asteroid type S,M,E,C,P,D"
 write(244,"(30(A12))") "#&
  S(0.5 km)", "S(1.0 km)","S(1.5 km)", "S(2.0 km)", "S(2.5 km)", &
 "M(0.5 km)", "M(1.0 km)","M(1.5 km)", "M(2.0 km)", "M(2.5 km)", &
 "E(0.5 km)", "E(1.0 km)","E(1.5 km)", "E(2.0 km)", "E(2.5 km)", &
 "C(0.5 km)", "C(1.0 km)","C(1.5 km)", "C(2.0 km)", "C(2.5 km)", &
 "P(0.5 km)", "P(1.0 km)","P(1.5 km)", "P(2.0 km)", "P(2.5 km)", &
 "D(0.5 km)", "D(1.0 km)","D(1.5 km)", "D(2.0 km)", "D(2.5 km)"

 do
 read(243,*,END=666) t,i,r,delta,alpha,lambda,beta

 alpha=alpha*180d0/pi


 ! TODO - Make the code neater
 ! FUNCTION IN: H, G1, G2, apha
 ! FUNCTION OUT: VIS1AU
 ! MAKE ANOTHER FUNCTION VIS1AU -> VISUAL - IN: VIS1AU, R, DELTA, OUT: VISUAL
 do k1=1,6
  do k2=1,5
   VIS1AU=10d0**(-0.4d0*H(k2,k1))*(G1(k1)*phi1(alpha)+G2(k1)*phi2(alpha)+(1d0-G1(k1)-G2(k1))*phi3(alpha))
   VIS1AU=log10(VIS1AU)/(-0.4d0)
   VISUAL=VIS1AU+5d0*log10(r*delta)

   VISLINE(nint((k1-1)*5d0+k2))=VISUAL
   enddo
  enddo
  
 write(244,"(30(F12.5))") VISLINE(1:30)
 enddo
 666 close(243)
 close(244)
 END
 
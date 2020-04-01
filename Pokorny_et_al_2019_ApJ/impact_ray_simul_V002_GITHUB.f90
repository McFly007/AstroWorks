! CHANGE LOG:
! VER 002 - Mar 31 2020 - by Petr Pokorny (petr.pokorny@nasa.gov / pokorny@cua.edu)
!
! Author:
! Petr Pokorny (petr.pokorny@nasa.gov)
! 
! LICENCE:
!                 ____________  
! This program is !!! FREE !!! to use, however, it would be very kind of you to
!                 ------------
! put a reference into your article (reference to Pokorny et al. (2020) - https://arxiv.org/abs/2003.12640 (soon published in the Astrophysical Journal)

! This program computes the meteoroid bombardment quantities onto a surface of a selected body (originally developed for the lunar poles)
! The code loads a file that describes the surface in triangles (i.e. obj file or a grid file that is sliced into triangles)
! and computes the incident angle of each meteoroid direction onto every surface triangle.
! A big portion of the code is dedicated to calculating whether each triangle is exposed to the meteoroid bombardment or it is "shadowed" by another surface triangle
! Once we have the incident angles for all triangles, we can add up the meteoroid quantities
! 
! Disclaimer: This code will be expanded, commented more thoroughly and made more efficient. Please do not hesitate to contact the author regarding 
! any comments or improvements
!
! Thank you for using this work


  module my_functions
    implicit none
    real*8, parameter :: pi=3.14159265358979323846d0
    real*8, parameter :: deg2rad=0.017453292520d0
    real*4,allocatable :: origtria(:,:) ! Array with Original Triangles, real*4 takes a lot less memory, the precision is enough
    real*4,allocatable :: rottria(:,:) ! Array with Rotated Triangles
    real*4,allocatable :: bbox(:,:) ! Array of Bounding Boxes
    real*4,allocatable :: triasort(:,:)   ! Array for Sorted Indices pointing to Rotated/Original Triangles
    real*4,allocatable :: lightC(:) ! Array of Incident angles for given triangles and the light/meteoroid ecliptic longitude/latitude
    real*4,allocatable :: lightM(:,:) ! Array of Meteoroid induced quantities for all triangles 
    real*8 rotmat(3,3),maxz ! Rotational matrix, and maxz of a triangle
    integer :: ndim,OFFS,ctr
    real*8 :: radmap2(180,90,3) ! Radiant distribution of meteoroid mass fluxes (col 1), meteoroids energy fluxes (col 2), and meteoroids ejecta production rate (col 3)
    real :: tnow=0.0,tprev=0.0
    parameter (ndim=1605,OFFS=5) ! Parameters for the original Moon maps
  contains

  recursive subroutine quicksort(a,col)
! This is a Quicksort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.)
! I modified it to arrays with multiple columns
    implicit none
    real*4 :: a(:,:) ! Array we want to sort
    integer :: col ! Column that will control the sorting
    real*4 x, t(col)
    integer :: first = 1, last
    integer i, j

    last = size(a, 1)
    x = a( (first+last) / 2 ,col)
    i = first
    j = last

    do
      do while (a(i,col) < x)
        i=i+1
      end do
      do while (x < a(j,col))
        j=j-1
      end do
      if (i >= j) exit
      t(:) = a(i,:);  a(i,:) = a(j,:);  a(j,:) = t(:)
      i=i+1
      j=j-1
    end do

    if (first < i - 1) call quicksort(a(first : i - 1,:),col)
    if (j + 1 < last)  call quicksort(a(j + 1 : last,:),col)
  end subroutine quicksort

  FUNCTION POINT_INSIDE_TRIANGLE(v1,v2,v3) result(vhit)
    ! Finds a centroid of a triangle in 3D
    implicit none
    real*8, intent(in) :: v1(3),v2(3),v3(3)
    real*8 :: s,t,vhit(3)
    integer :: cntr
    cntr=0
    do while(cntr.lt.1d0)
      s=0.33333d0
      t=0.33333d0
      s=1d0/3d0
      t=1d0/3d0
      if((s+t).le.1) then 
        vhit=v1+s*(v2-v1)+t*(v3-v1)
        cntr=cntr+1
      endif
    enddo
  END FUNCTION POINT_INSIDE_TRIANGLE
 
  FUNCTION POINT2D(v1,v2,v3) result(vhit)
    ! Finds a centroid of a triangle in 2D
    implicit none
    real*4, intent(in) :: v1(2),v2(2),v3(2)
    real*4 :: s,t,ran3,vhit(2)
    integer :: cntr
    cntr=0
    do while(cntr.lt.1d0)
      s=0.33333d0
      t=0.33333d0
      s=1d0/3d0
      t=1d0/3d0
      if((s+t).le.1) then 
        vhit=v1+s*(v2-v1)+t*(v3-v1)
        cntr=cntr+1
      endif
    enddo
  END FUNCTION POINT2D
 
  PURE FUNCTION scalar(a, b)
  ! Scalar/dot product of a . b in 3D
    REAL*8 :: scalar
    REAL*8, DIMENSION(3), INTENT(IN) :: a, b
    scalar = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
  END FUNCTION scalar
 
  FUNCTION scalar2d(a, b)
  ! Scalar/dot product of a . b in 2D
    REAL*8 :: scalar2d
    REAL*8, DIMENSION(2), INTENT(IN) :: a, b
    scalar2d = a(1)*b(1)+a(2)*b(2)
  END FUNCTION scalar2d
 
  PURE FUNCTION cross(a, b)
  ! Cross product of a x b in 3D  
    REAL*8, DIMENSION(3) :: cross
    REAL*8, DIMENSION(3), INTENT(IN) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION cross
  
  FUNCTION CLEANSPACE(flnm) result(flnm2)
    ! Cleans up strings from blank spaces
    implicit none
    character*50 :: flnm,flnm2
    integer :: i
    flnm2=""
    do i=1,len(flnm)
      if (flnm(i:i).ne.' ')flnm2=trim(flnm2)//trim(flnm(i:i))   
    end do
  END FUNCTION CLEANSPACE

  PURE FUNCTION INSIDE_TRIA(v1,v2,v3,rayorig,rayvec,rimin) result(ri)
  ! Checks if a ray originating in rayorig, defined by a vector rayvec is passing through triangle (v1,v2,v3)
  ! Result: if yes, then "ri" is the distance from rayorig to a point inside the triangle (v1,v2,v3)\
  !         if not, then "ri = -9d9", so we can check if the result is negative afterwards 
  !         (distance should be always positive if the ray hits the triangle)
  ! I'm using this http://geomalgorithms.com/a06-_intersect-2.html
    implicit none
    real*8,intent(in) :: rayorig(3),rayvec(3),rimin
    real*8,intent(in) :: v1(3),v2(3),v3(3)
    real*8, dimension(3) :: v1v2,v1v3,vnorm,rivec1,rivec2
    real*8, dimension(3) :: vecw,vecu,vecv
    real*8 :: ri,si,ti,suv,svw,svv,suw,suu,denom
    v1v2=v1-v2
    v1v3=v1-v3
    vnorm=cross(v1v2,v1v3)
    rivec1=v1-rayorig
    rivec2=-rayvec
    ri=scalar(vnorm,rivec1)/scalar(vnorm,rivec2)
    vecw=rayorig-ri*rayvec-v1
    vecu=v3-v1
    vecv=v2-v1

    suv=scalar(vecu,vecv)
    svw=scalar(vecw,vecv)
    svv=scalar(vecv,vecv)
    suw=scalar(vecw,vecu)
    suu=scalar(vecu,vecu)

    si=suv*svw-svv*suw
    ti=suv*suw-suu*svw
    denom=1d0/(suv*suv-suu*svv)

    si=si*denom
    ti=ti*denom
    if((ri.lt.(rimin-1d-5)).and.(ri.gt.0d0).and.(si.ge.0d0)&
      .and.(ti.ge.0d0).and.((si+ti).le.1d0)) then
      ri=ri
    else
      ri=-9d9  
    endif
  END FUNCTION INSIDE_TRIA

  PURE FUNCTION GET_RAYVEC(lon,lat) result(rayvec)
  ! Calculates the rayvec based on the raylongitude and latitude in radians
    implicit none
    real*8,intent(IN) :: lon,lat
    real*8,dimension(3) :: rayvec
    rayvec(1)=cos(lon)*cos(lat)
    rayvec(2)=sin(lon)*cos(lat)
    rayvec(3)=sin(lat)
  END FUNCTION GET_RAYVEC

  PURE FUNCTION GET_ROTMAT(lon,lat) result(rotmat)
  ! Calculates the rotational matrix based on the ray longitude and latitude
    implicit none
    real*8,intent(IN) :: lon,lat
    real*8, dimension (3) :: rotz1,rotz2,rotz3
    real*8, dimension (3) :: roty1,roty2,roty3
    real*8, dimension (3,3) :: rotz,roty,rotmat

    rotz1=[dcos(lon),-dsin(lon),0d0]
    rotz2=[dsin(lon), dcos(lon),0d0]
    rotz3=[0d0,0d0,1d0]

    roty1=[dcos(0.5d0*pi-lat),0d0,dsin(0.5d0*pi-lat)]
    roty2=[0d0, 1d0, 0d0]
    roty3=[-dsin(0.5d0*pi-lat),0d0,dcos(0.5d0*pi-lat)]
    rotz(:,1)=rotz1
    rotz(:,2)=rotz2
    rotz(:,3)=rotz3

    roty(:,1)=roty1
    roty(:,2)=roty2
    roty(:,3)=roty3

    rotmat=matmul(roty,rotz)
  END FUNCTION GET_ROTMAT   

  SUBROUTINE READ_GRIDFILE(flnm)
    ! Reads a binary grid file used to calculate results in Pokorny et al. (2020)
    ! Fills "origtria" and allocates "rottria" and "bbox" and calculates the total number of triangles "ctr"
    implicit none
    character*50 flnm
    real*4 prev
    integer n1,n2,i,j,i1,j1
    real*4,allocatable :: xy(:,:,:)
    real*4 reado(3)

    allocate(xy(ndim,ndim,3))
    xy=0d0
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old",form="unformatted")
    write(*,*) "Reading the following grid file:",flnm
    i=1
    j=1
    do 
      read(44,end=666) reado(1),reado(2),reado(3)
      if(reado(2).lt.prev) then
        if(n2.lt.j) n2=j
        i=i+1
        j=2
      endif
      xy(i,j,:)=reado(:)
      prev=reado(2)
      j=j+1
    enddo

    666 close(44)
    n1=i-2
    n2=j-2

    ctr=0
    allocate(origtria((n1*n2)*2,10))
    do i1=1+OFFS,n1-OFFS-1
      do j1=1+OFFS,n2-OFFS-1
        ctr=ctr+1
        origtria(ctr,1:9)=(/xy(i1,j1,:),xy(i1+1,j1,:),xy(i1,j1+1,:)/)
      enddo
    enddo
    write(*,*) "FIRST SET OF TRIANGLES:",ctr
    do i1=n1-OFFS,2+OFFS,-1
      do j1=n2-OFFS,2+OFFS,-1
        ctr=ctr+1
        origtria(ctr,1:9)=(/xy(i1,j1,:),xy(i1-1,j1,:),xy(i1,j1-1,:)/)
      enddo
    enddo
    write(*,*) "SECOND SET OF TRIANGLES",ctr

    deallocate(xy)
    allocate(rottria(ctr,10))
    allocate(bbox(ctr,4))
  END SUBROUTINE READ_GRIDFILE

  SUBROUTINE READ_OBJFILE(flnm)
    ! This subroutine reads obj files that are defined by vectors "v" and faces "f"
    ! Fills "origtria" and allocates "rottria" and "bbox" and calculates the total number of triangles "ctr"
    implicit none
    character*50 flnm
    character*1 idstr
    integer veccount,facecount
    real*4 reado(3)
    integer intdo(3),i
    real*4,allocatable ::  vecarr(:,:)
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old")
    write(*,*) flnm
    veccount=0
    facecount=0

    do 
      read(44,*,end=666) idstr,reado(1),reado(2),reado(3)
      if(idstr.eq."v") then
        veccount=veccount+1
      elseif(idstr.eq."f") then
        facecount=facecount+1
      endif  
    enddo
    666 close(44)   
    open(44,file=trim(flnm),status="old")
    allocate(origtria(facecount,10))
    allocate(vecarr(veccount,3))
    origtria=0d0
    do i=1,veccount+facecount
      read(44,*) idstr,reado(1),reado(2),reado(3)
      ! write(*,*) idstr,reado(1),reado(2),reado(3)
      if(idstr.eq."v") then
        vecarr(i,:)=reado
      elseif(idstr.eq."f") then
        intdo(1)=nint(reado(1))
        intdo(2)=nint(reado(2))
        intdo(3)=nint(reado(3))
        ! facecount=facecount+1
        origtria(i-veccount,1:9)=(/vecarr(intdo(1),:),vecarr(intdo(2),:),vecarr(intdo(3),:)/)
      endif  
    enddo
    ! 666 close(44)   
    ctr=facecount
    allocate(bbox(ctr,4))
    allocate(rottria(ctr,10))
    deallocate(vecarr)
    write(*,*) "READ OBJFILE DONE",ctr
  END SUBROUTINE READ_OBJFILE

  SUBROUTINE READ_BINOBJFILE(flnm)
    ! This subroutine reads precalculated binary obj files (which is faster than reading ascii files)
    implicit none
    character*50 flnm
    character*1 idstr
    integer veccount,facecount
    real*8 reado(3)
    integer intdo(3)
    real*4,allocatable :: tmparr(:,:)
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old",form="unformatted")
    write(*,*) flnm

    read(44) veccount,facecount
    ctr=facecount
    allocate(origtria(ctr,10))
    allocate(bbox(ctr,4))
    allocate(rottria(ctr,10))
    allocate(tmparr(ctr,10))
    read(44) tmparr

    666 close(44)   
    origtria=tmparr
    deallocate(tmparr)
    ctr=facecount
    write(*,*) "READ BINOBJFILE DONE",ctr
  END SUBROUTINE READ_BINOBJFILE

  SUBROUTINE LOAD_METEOROIDS()
    ! This subroutine loads meteoroid files and populates "radmap2"
    implicit none
    real*8 lon,lat,vel,flux,tons(4)
    integer ilon,ilat,ivel,i
    open(21, file="./AST", status="OLD")
    open(22, file="./JFC", status="OLD")
    open(23, file="./HTC", status="OLD")
    open(24, file="./OCC", status="OLD")
    write(*,*) "How many tons of AST,JFC,HTC,OCC (in tons, comma separated)?"
    ! read(*,*) tons
    tons=(/3.7,34.3,2.82,2.18/) ! From Pokorny et al. (2019) - Meteoroids at the Moon in tons/day

    do i=1,4
      do
        read(20+i,*,end=555) lon,lat,vel,flux 
        ilon=nint((lon+181d0)/2d0)
        ilat=nint((lat+91d0)/2d0)
        ivel=nint((vel+1d0)/2d0)
        radmap2(ilon,ilat,1)=radmap2(ilon,ilat,1)+flux*tons(i)
        radmap2(ilon,ilat,2)=radmap2(ilon,ilat,2)+flux*tons(i)*vel**2*0.5 !ENERGY 1/2*m*v^2
        radmap2(ilon,ilat,3)=radmap2(ilon,ilat,3)+flux*tons(i)*vel**2.46*7.358 ! WITH K&G2001 factor 
      enddo
      555 close(20+i)
    enddo

  END SUBROUTINE LOAD_METEOROIDS
 
  SUBROUTINE TIMER()
    ! This little routine tracks time for us
    CALL CPU_TIME(tnow)
    write(*,*) "TOTAL RUNTIME:",tnow,"SEGMENT RUNTIME",tnow-tprev
    tprev=tnow
  END SUBROUTINE TIMER

  SUBROUTINE OPEN_OUTPUT(lonbeg,lonend,latbeg,latend)
    ! Creates an output binary file using the initial parameters
    implicit none
    real*8 lonbeg,lonend,lonstep,latbeg,latend
    character*50 :: flnm2=""

    write(flnm2,"(A6,I4,A1,I4,A1,I4,A1,I4)") "IMPACTS_",nint(lonbeg),"_",nint(lonend),"_",nint(latbeg),"_",nint(latend)
    flnm2=CLEANSPACE(flnm2)
    write(*,*) "The output file is:",flnm2
    open(66,file=trim(flnm2),status="unknown",form="unformatted")
  END SUBROUTINE OPEN_OUTPUT

end module my_functions
 
! End of modules, now the real program starts here !
 

 
PROGRAM MAPME

  use my_functions
  implicit none
  include 'vars.inc'
  
  real*8 v1(3),v2(3),v3(3),vnorm(3),rnd,s,t,ran3,surfvec(3),surf
  real*8 rotvhit(3),cosdec
  real*8 lon,lat,rayvec(3),rayorig(3),ri,rimin,vecw(3)
  real*8 si,ti,vecu(3),vecv(3),denom,radius,vhit(3),prev
  real*8 rotz1(3),rotz2(3),rotz3(3),v1v2(3),v1v3(3),rivec1(3),rivec2(3)
  real*8 rotx1(3),rotx2(3),rotx3(3)
  real*8 roty1(3),roty2(3),roty3(3)
  real*8 vrot1(3),vrot2(3),vrot3(3),suu,suv,svv,suw,svw
  real*8 rotz(3,3), rotx(3,3),roty(3,3)
  integer i,j,n,k,i1,j1,cntrmax,kkk,n1,n2
  integer ilon,ilat
  logical :: insider
  real*8 :: rlim=1d0
  character*50 :: flnm="", flnm2=""

  real*8 lonbeg,lonend,lonstep,latbeg,latend,latstep
  integer kk1,kk2

  ! 
  write(*,*) "Enter the min longitude (deg), max longitude (deg), longitude step (deg),&
   min latitude (deg), max latitude (deg), latitude step (deg), triangle file name"
  read(*,*) lonbeg,lonend,lonstep,latbeg,latend,latstep,flnm

  ! call READ_BINOBJFILE(flnm)
  call READ_GRIDFILE(flnm)
  call OPEN_OUTPUT(lonbeg,lonend,latbeg,latend)
  call LOAD_METEOROIDS() !!! FILLS RAD2MAP with numbers

  allocate(triasort(ctr,2))
  allocate(lightM(ctr,3))
  allocate(LightC(ctr))

  lightM=0d0

  do kk1=nint(lonbeg),nint(lonend),nint(lonstep) ! We iterate over the entire longitude range 
    do kk2=nint(latbeg),nint(latend),nint(latstep) ! We iterate over the entire latitude range - These two loops are independent and can be split into multiple serial runs
      lightC=0d0 
      write(*,*) "STARTED, LON:",kk1,"LAT:",kk2
      lon=kk1
      lat=kk2   
      ilon=nint((lon+181d0)/2d0)
      ilat=nint((lat+91d0)/2d0)

      lon=lon*deg2rad
      lat=lat*deg2rad

      rotmat=GET_ROTMAT(lon,lat)
      rayvec=GET_RAYVEC(lon,lat)

      do i1=1,ctr
        rottria(i1,1:9)=(/matmul(rotmat,origtria(i1,1:3)),matmul(rotmat,origtria(i1,4:6)),matmul(rotmat,origtria(i1,7:9))/)  ! rotate the array of triangles so it the ray is coming from z-axis direction and the problem is actually 2D
        maxz=maxval(rottria(i1,3:9:3)) ! Get the max z values for each triangle 
        triasort(i1,1)=i1
        triasort(i1,2)=maxz
        bbox(i1,1:2)=POINT2D(rottria(i1,1:2),rottria(i1,4:5),rottria(i1,7:8)) ! Prepare the bounding boxes 
        bbox(i1,3)= maxval((/(abs(bbox(i1,1:2)-rottria(i1,1:2))),&
          (abs(bbox(i1,1:2)-rottria(i1,4:5))),&
          (abs(bbox(i1,1:2)-rottria(i1,7:8)))/))
        bbox(i1,4)=maxz
      enddo

      call quicksort(triasort(1:ctr,:),2) ! Sort triangles by their rotated z-axis value
      call quicksort(bbox(1:ctr,:),4) ! Sort bounding boxes

      ! Create the bounding box grid for faster lookup (this will be replaced by voxels and AABB trees)
      r_boxminX=minval(bbox(1:ctr,1))-maxval(bbox(1:ctr,3))
      r_boxmaxX=maxval(bbox(1:ctr,1))+maxval(bbox(1:ctr,3))
      r_boxminY=minval(bbox(1:ctr,2))-maxval(bbox(1:ctr,3))
      r_boxmaxY=maxval(bbox(1:ctr,2))+maxval(bbox(1:ctr,3))
      r_dx=(r_boxmaxX-r_boxminX)/(r_num*1d0)
      r_dy=(r_boxmaxY-r_boxminY)/(r_num*1d0)

      r_pivots=0
      r_numbs=0

      do i=1,ctr
        r_id(1)=int((bbox(i,1)-bbox(i,3)-r_boxminX)/r_dx)+1
        r_id(2)=int((bbox(i,1)+bbox(i,3)-r_boxminX)/r_dx)+1
        r_id(3)=int((bbox(i,2)-bbox(i,3)-r_boxminY)/r_dY)+1
        r_id(4)=int((bbox(i,2)+bbox(i,3)-r_boxminY)/r_dY)+1
        do l1=1,4
          if(r_id(l1).le.0) r_id(l1)=1
          if(r_id(l1).gt.r_num) r_id(l1)=r_num
        enddo



        do l1=r_id(1),r_id(2)
          do l2=r_id(3),r_id(4)
            r_pivots((l1-1)*r_num+l2)=r_pivots((l1-1)*r_num+l2)+1
          enddo
        enddo


      enddo
      ! Here we create pivots for faster access and shorter processing
      allocate(r_tria(sum(r_pivots),11))
      r_tria=0d0

      do l1=1,r_num*r_num
        do l2=1,l1-1
          r_numbs(l1)=r_numbs(l1)+r_pivots(l2)
        enddo
        r_numbs(l1)=r_numbs(l1)+1
      enddo

      do i=1,ctr
        r_id(1)=int((bbox(i,1)-bbox(i,3)-r_boxminX)/r_dx)+1
        r_id(2)=int((bbox(i,1)+bbox(i,3)-r_boxminX)/r_dx)+1
        r_id(3)=int((bbox(i,2)-bbox(i,3)-r_boxminY)/r_dY)+1
        r_id(4)=int((bbox(i,2)+bbox(i,3)-r_boxminY)/r_dY)+1
        do l1=1,4
          if(r_id(l1).le.0) r_id(l1)=1
          if(r_id(l1).gt.r_num) r_id(l1)=r_num
        enddo

        do l1=r_id(1),r_id(2)
          do l2=r_id(3),r_id(4)
            r_index=r_numbs((l1-1)*r_num+l2)
            r_numbs((l1-1)*r_num+l2)=r_numbs((l1-1)*r_num+l2)+1
            r_tria(r_index,1:10)=rottria(i,:)
            r_tria(r_index,11)=i
          enddo
        enddo


      enddo

      ! Now we start determining if all triangles are shadowed and what is the incident angle
      radius=100000d0 ! rays start somewhere far

      do i1=1,ctr
        v1=origtria(nint(triasort(i1,1)),1:3) ! We are going through all triangles based on their rotated z-axis distance
        v2=origtria(nint(triasort(i1,1)),4:6)
        v3=origtria(nint(triasort(i1,1)),7:9)
        v1v2=v2-v1
        v1v3=v3-v1
        surfvec=cross(v1v2,v1v3)
        if(SCALAR(surfvec,v1).lt.0) surfvec=-surfvec
        if(SCALAR(rayvec,surfvec).lt.0) THEN
          cycle
        ENDIF
        cosdec=SCALAR(rayvec,surfvec)/dsqrt(scalar(surfvec,surfvec))  ! Rayvec is normalized
        !surf=scalar(surfvec,surfvec)*0.5d0 ! Calcuates the surface area

        vhit=POINT_INSIDE_TRIANGLE(v1,v2,v3)
        v1=rottria(nint(triasort(i1,1)),1:3)
        v2=rottria(nint(triasort(i1,1)),4:6)
        v3=rottria(nint(triasort(i1,1)),7:9)
        rotvhit=POINT_INSIDE_TRIANGLE(v1,v2,v3)

        rayorig=vhit+radius*rayvec
        rimin=radius

        r_id(1)=int((rotvhit(1)-r_boxminX)/r_dx)+1
        r_id(2)=int((rotvhit(2)-r_boxminY)/r_dY)+1

        ! This loop goes from all triangles that can shadow a selected triangle - we calculate if the shadowing occurs using INSIDE_TRIA function
        do l1=r_numbs( (r_id(1)-1)*r_num+r_id(2))-r_pivots( (r_id(1)-1)*r_num+r_id(2)),&
          r_numbs( (r_id(1)-1)*r_num+r_id(2))
          i=nint(r_tria(l1,11))

          v1=origtria(nint(triasort(i,1)),1:3)
          v2=origtria(nint(triasort(i,1)),4:6)  
          v3=origtria(nint(triasort(i,1)),7:9)
          ri=INSIDE_TRIA(v1,v2,v3,rayorig,rayvec,rimin)
          if(ri.gt.0d0) then
            rimin=ri
            goto 881 ! If we have positive value then the triangle is shadowed => break the loop and continue
          endif
        enddo

        881 continue
        if(rimin.gt.(radius-1e-5)) then ! If nothing is shadowed rimin = radius and then we put cosdec into LightC
          lightC(nint(triasort(i1,1)))=cosdec
        endif
      enddo


      ! Here we add the meteoroid quantities - LightC has the cosine of the indicent angle (a.k.a cosdec)
      do ilon=1,180 ! Here we assume that Moon rotates around the z-axis (i.e. the sky rotates along the z-axis) and we assume that the moon rotates exactly once (i.e. we average over one rotation of the Moon)
        do i=1,ctr
          lightM(i,1)=lightM(i,1)+lightC(i)*radmap2(ilon,ilat,1)
          lightM(i,2)=lightM(i,2)+lightC(i)*radmap2(ilon,ilat,2)
          lightM(i,3)=lightM(i,3)+lightC(i)**3d0*radmap2(ilon,ilat,3)
        enddo
      enddo
      call TimeR()

      deallocate(r_tria)
    enddo
  enddo
  lightM=lightM/180d0 ! We added up 180 different sky rotations, so we divide by that number to get an average
  lightM=lightM/(637800000d0**2d0*pi)  ! This transforms the results to cm^2 (the original meteoroid files are in kg per day per cross-section of Earth)
  lightM=lightM*1d3/(86400d0) ! Transform to g/cm^2/s
  do i=1,ctr
    write(66) real(origtria(i,1:3)),real(lightM(i,1:3))
  enddo
  1235 format (2(1x,f9.5),1(1x,f11.5))
end
 
 
  function ran3(idum)
! Random number generator - (real*8 version based on Numerical Recipes).
  implicit real*8(a-h,o-z)
  parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.d-9)
  dimension ma(55)
  data iff /0/
  save

  if((idum.ge.0).and.(iff.ne.0)) goto 20
    iff=1
    mj=mseed-iabs(idum)
    mj=mod(mj,mbig)
    ma(55)=mj
    mk=1
    do 11 i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if(mk.lt.mz)mk=mk+mbig
      mj=ma(ii)
11  continue
    do 13 k=1,4
      do 12 i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if(ma(i).lt.mz)ma(i)=ma(i)+mbig
12    continue
13  continue
    inext=0
    inextp=31
    idum=1
20  continue
  inext=inext+1
  if(inext.eq.56)inext=1
  inextp=inextp+1
  if(inextp.eq.56)inextp=1
  mj=ma(inext)-ma(inextp)
  if(mj.lt.mz)mj=mj+mbig
  ma(inext)=mj
  ran3=mj*fac
  return
  end
      
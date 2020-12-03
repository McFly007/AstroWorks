! Functions and Subroutines for Meteoroids ray-tracing code for OBJ files using FASTBVH and TINYOBJ
! Written by Petr Pokorny (petr.pokorny@nasa.gov / pokorny@cua.edu) 
!
! Disclaimer: This code is not final (but works) and will be properly commented soon.
! V03 - Dec 1st, 2020
module my_functions
    use omp_lib
    implicit none
    real*4, parameter :: pi=3.14159265358979323846e0
    real*4, parameter :: deg2rad=0.017453292520e0
    real*4, parameter :: rad2deg=57.295779513e0
    real*4,allocatable :: origtria(:,:) ! Array with Original Triangles, real*4 takes a lot less memory, the precision is enough
    real*4,allocatable :: rottria(:,:) ! Array with Rotated Triangles
    real*4,allocatable :: bbox(:,:) ! Array of Bounding Boxes
    real*4,allocatable :: triasort(:,:)   ! Array for Sorted Indices pointing to Rotated/Original Triangles
    real*4,allocatable :: lightC(:) ! Array of Incident angles for given triangles and the light/meteoroid ecliptic longitude/latitude
    real*4,allocatable :: lightM(:,:) ! Array of Meteoroid induced quantities for all triangles
    real*4,allocatable :: lightMPar(:,:,:) ! Array of Meteoroid induced quantities for all triangles 
    real*4 rotmat(3,3),maxz ! Rotational matrix, and maxz of a triangle
    real*8 orbel_vec(5)
    integer :: ndim,OFFS,ctr
    integer :: NFILES ! number of meteoroid files that are processed
    integer :: comp_ctr ! Number of faces where the code does calculations (smaller for specialized area calculations)
    real*4 :: radmap2(180,90,3) ! Radiant distribution of meteoroid mass fluxes (col 1), meteoroids energy fluxes (col 2), and meteoroids ejecta production rate (col 3)
    real*4,allocatable :: radmap2par(:,:,:,:),radmapAVG(:,:,:)
    real, allocatable :: ARvhit(:,:),ARsurfvec(:,:),elon(:)
    real*4, allocatable :: TAA(:)
    real, allocatable :: impactC(:),radmapTOTAL(:,:,:,:)
    real :: tnow=0.0,tprev=0.0
    real*4 :: tilt ! The axial tilt
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
    real*4, intent(in) :: v1(3),v2(3),v3(3)
    real*4 :: s,t,vhit(3)
    integer :: cntr
    cntr=0
    do while(cntr.lt.1e0)
      s=0.33333e0
      t=0.33333e0
      s=1e0/3e0
      t=1e0/3e0
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
    do while(cntr.lt.1e0)
      s=0.33333e0
      t=0.33333e0
      s=1e0/3e0
      t=1e0/3e0
      if((s+t).le.1) then 
        vhit=v1+s*(v2-v1)+t*(v3-v1)
        cntr=cntr+1
      endif
    enddo
  END FUNCTION POINT2D
 
  PURE FUNCTION scalar(a, b)
  use omp_lib
  ! Scalar/dot product of a . b in 3D
    real*4 :: scalar
    real*4, DIMENSION(3), INTENT(IN) :: a, b
    scalar = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
  END FUNCTION scalar
 
  FUNCTION scalar2d(a, b)
  ! Scalar/dot product of a . b in 2D
    real*4 :: scalar2d
    real*4, DIMENSION(2), INTENT(IN) :: a, b
    scalar2d = a(1)*b(1)+a(2)*b(2)
  END FUNCTION scalar2d
 
  PURE FUNCTION cross(a, b)
  ! Cross product of a x b in 3D  
    real*4, DIMENSION(3) :: cross
    real*4, DIMENSION(3), INTENT(IN) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION cross
  
  FUNCTION CLEANSPACE(flnm) result(flnm2)
    ! Cleans up strings from blank spaces
    implicit none
    character*150 :: flnm,flnm2
    integer :: i
    flnm2=""
    do i=1,len(flnm)
      if (flnm(i:i).ne.' ')flnm2=trim(flnm2)//trim(flnm(i:i))   
    end do
  END FUNCTION CLEANSPACE

  PURE FUNCTION GET_ROTMATXYZ(angle,which) result(rotmatxyz)
  ! Calculates the rotational matrix based on the ray longitude and latitude
    implicit none
    real*4,intent(IN) :: angle
    real*4, dimension (3) :: rotz1,rotz2,rotz3
    real*4, dimension (3) :: roty1,roty2,roty3
    real*4, dimension (3) :: rotx1,rotx2,rotx3
    real*4, dimension (3,3) :: rotz,roty,rotmatxyz,rotx
    character*1,intent(IN):: which

    real*4 :: pi
    pi =acos(-1e0)
    rotz1=[cos(angle),-sin(angle),0e0]
    rotz2=[sin(angle), cos(angle),0e0]
    rotz3=[0e0,0e0,1e0]

    roty1=[cos(angle),0e0,sin(angle)]
    roty2=[0e0, 1e0, 0e0]
    roty3=[-sin(angle),0e0,cos(angle)]

    rotx1=[1e0, 0e0, 0e0]
    rotx2=[0e0,cos(angle),-sin(angle)]
    rotx3=[0e0,sin(angle),cos(angle)]

    rotz(1,:)=rotz1
    rotz(2,:)=rotz2
    rotz(3,:)=rotz3

    roty(1,:)=roty1
    roty(2,:)=roty2
    roty(3,:)=roty3

    rotx(1,:)=rotx1
    rotx(2,:)=rotx2
    rotx(3,:)=rotx3

    if(which.eq."X") rotmatxyz=rotx
    if(which.eq."Y") rotmatxyz=roty
    if(which.eq."Z") rotmatxyz=rotz

  END FUNCTION GET_ROTMATXYZ   

    FUNCTION GET_SPHERE(lon_in,lat_in,tilt,iphi)
  ! lon,lat,tile are in radians
  ! iphi is between 0 and 358
  ! Tilt is the axial tilt
  ! Ceres has - α = 291.42751 ∘ and δ = 66.76043 ∘ - spin axis in RA and Dec
  ! Which is Lon = 11.186218 Lat = 81.550375, ceres inclination is 10.6 degrees
! 11.18702682 degrees
! Output Latitude or DEC (degrees)  81.55035527 degrees
    real*4,intent(IN) :: lon_in,lat_in,tilt
    integer,intent(IN) :: iphi
    real*4 rho,r,phi,rand,lat
    real*4 ovec(3),fvec(3),rmat(3,3),v1(3),v2(3),plane_vec(3),center(3),orig_center(3),rad
    integer i,ilon,ilat
    integer :: GET_SPHERE(2)
    lat=lat_in
    ! This is just to secure that if we have tilt equal to lat, we do not get NaNs 
    if(abs(lat-tilt).lt.1e-2) then
        call random_number(rand)
        lat=lat+1e-2*(-1.0**int(rand*2.0))
    endif    

    rho=sqrt(1e0-cos(lat-tilt)**2.0)

    rmat=GET_ROTMATXYZ(tilt,"Y")
    orig_center=[0e0,0e0,rho]
    plane_vec=matmul(rmat,orig_center)
    ovec=[cos(lon_in)*cos(lat),sin(lon_in)*cos(lat),sin(lat)]
    center=scalar(ovec,plane_vec)/(scalar(plane_vec,plane_vec))*plane_vec

   
    v1=ovec-center
    rad=sqrt(scalar(v1,v1))
    v2=cross(plane_vec,v1)
    ! Here we calculate the celestial circle - ilon and ilat are integer values from 1-180 and 1-90
    ! do i=0,358,2
        phi=real(iphi)*deg2rad
        fvec=center+rad*cos(phi)*v1/sqrt(scalar(v1,v1))+rad*sin(phi)*v2/sqrt(scalar(v2,v2))
        ilon=int((atan2(fvec(2),fvec(1))/deg2rad+180e0)*0.5)+1
        ilat=int((asin(fvec(3))/deg2rad+90e0)*0.5)+1
        ! write(*,*) fvec,ilon,ilat,lon_in,lat_in
    ! enddo
        GET_SPHERE=[ilon,ilat]


 END FUNCTION GET_SPHERE

    FUNCTION GET_SPHERE2(lon_in,lat_in,tilt,elon_in,iphi)
  ! lon,lat,tilt,elon are in radians
  ! iphi is between 0 and 358
  ! Tilt is the axial tilt
  ! Elon is the angle between the apex of the motion of the lon of the spin axis
  ! Ceres has - α = 291.42751 ∘ and δ = 66.76043 ∘ - spin axis in RA and Dec
  ! Which is Lon = 11.186218 Lat = 81.550375, ceres inclination is 10.6 degrees
    real*4,intent(IN) :: lon_in,lat_in,tilt,elon_in
    integer,intent(IN) :: iphi
    real*4 rho,r,phi,rand,lat,fvec2(3)
    real*4 ovec(3),fvec(3),rmat(3,3),v1(3),v2(3),plane_vec(3),center(3),orig_center(3),rad
    integer i,ilon,ilat
    integer :: GET_SPHERE2(2)
    lat=lat_in
    ! This is just to secure that if we have tilt equal to lat, we do not get NaNs 
    if(abs(lat-tilt).lt.1e-2) then
        call random_number(rand)
        lat=lat+1e-2*(-1.0**int(rand*2.0))
    endif    

    ! if(abs(lat-elon_in).lt.1e-2) then
    !     call random_number(rand)
    !     lat=lat+1e-2*(-1.0**int(rand*2.0))
    ! endif    

    rho=sqrt(1e0-cos(lat-tilt)**2.0)
    ! rmat=GET_ROTMATXYZ(tilt,"Y")
    rmat=matmul(GET_ROTMATXYZ(tilt,"Y"),GET_ROTMATXYZ(elon_in,"Z"))
    ! rmat=matmul(GET_ROTMATXYZ(elon_in,"Z"),GET_ROTMATXYZ(tilt,"Y"))
    orig_center=[0e0,0e0,rho]
    plane_vec=matmul(rmat,orig_center)
    ovec=[cos(lon_in)*cos(lat),sin(lon_in)*cos(lat),sin(lat)]
    center=scalar(ovec,plane_vec)/(scalar(plane_vec,plane_vec))*plane_vec

    
    v1=ovec-center
    rad=sqrt(scalar(v1,v1))
    v2=cross(plane_vec,v1)
    ! write(*,*) rho,v1,v2,rad,elon_in/deg2rad,tilt/deg2rad
    ! Here we calculate the celestial circle - ilon and ilat are integer values from 1-180 and 1-90
    ! do i=0,358,2
        phi=real(iphi)*deg2rad
        fvec=center+rad*cos(phi)*v1/sqrt(scalar(v1,v1))+rad*sin(phi)*v2/sqrt(scalar(v2,v2))
        ! fvec2=matmul(GET_ROTMATXYZ(elon_in,"Z"),fvec)
        fvec2=fvec/sqrt(scalar(fvec,fvec))
        ilon=int((atan2(fvec2(2),fvec2(1))/deg2rad+180e0)*0.5)+1
        ilat=int((asin(fvec2(3))/deg2rad+90e0)*0.5)+1
         ! write(*,*) fvec,ilon,ilat,lon_in,lat_in
    ! enddo
        GET_SPHERE2=[ilon,ilat]


 END FUNCTION GET_SPHERE2
  PURE FUNCTION INSIDE_TRIA(v1,v2,v3,rayorig,rayvec,rimin) result(ri)
  ! Checks if a ray originating in rayorig, defined by a vector rayvec is passing through triangle (v1,v2,v3)
  ! Result: if yes, then "ri" is the distance from rayorig to a point inside the triangle (v1,v2,v3)\
  !         if not, then "ri = -9d9", so we can check if the result is negative afterwards 
  !         (distance should be always positive if the ray hits the triangle)
  ! I'm using this http://geomalgorithms.com/a06-_intersect-2.html
    implicit none
    real*4,intent(in) :: rayorig(3),rayvec(3),rimin
    real*4,intent(in) :: v1(3),v2(3),v3(3)
    real*4, dimension(3) :: v1v2,v1v3,vnorm,rivec1,rivec2
    real*4, dimension(3) :: vecw,vecu,vecv
    real*4 :: ri,si,ti,suv,svw,svv,suw,suu,denom
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
    denom=1e0/(suv*suv-suu*svv)

    si=si*denom
    ti=ti*denom
    if((ri.lt.(rimin-1e-5)).and.(ri.gt.0e0).and.(si.ge.0e0)&
      .and.(ti.ge.0e0).and.((si+ti).le.1e0)) then
      ri=ri
    else
      ri=-9d9  
    endif
  END FUNCTION INSIDE_TRIA

  PURE FUNCTION GET_RAYVEC(lon,lat) result(rayvec)
  ! Calculates the rayvec based on the raylongitude and latitude in radians
    implicit none
    real*4,intent(IN) :: lon,lat
    real*4,dimension(3) :: rayvec
    rayvec(1)=cos(lon)*cos(lat)
    rayvec(2)=sin(lon)*cos(lat)
    rayvec(3)=sin(lat)
  END FUNCTION GET_RAYVEC

  SUBROUTINE READ_GRIDFILE(flnm)
    ! Reads a binary grid file used to calculate results in Pokorny et al. (2020)
    ! Fills "origtria" and allocates "rottria" and "bbox" and calculates the total number of triangles "ctr"
    implicit none
    character*150 flnm
    real*4 prev
    integer n1,n2,i,j,i1,j1
    real*4,allocatable :: xy(:,:,:)
    real*4 reado(3)

    allocate(xy(ndim,ndim,3))
    xy=0e0
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
    character*150 flnm
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
    origtria=0e0
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
    character*150 flnm
    character*1 idstr
    integer veccount,facecount
    real*4 reado(3)
    integer intdo(3)
    real*4,allocatable :: tmparr(:,:)
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old",form="unformatted")
    write(*,*) flnm

    read(44) veccount,facecount
    ctr=facecount
    allocate(origtria(ctr,10))
    ! allocate(bbox(ctr,4))
    ! allocate(rottria(ctr,10))
    allocate(tmparr(ctr,10))
    read(44) tmparr

    666 close(44)   
    origtria=tmparr
    deallocate(tmparr)
    ctr=facecount
    write(*,*) "READ BINOBJFILE DONE",ctr
  END SUBROUTINE READ_BINOBJFILE

    SUBROUTINE READ_BINOBJFILE_COMBE(flnm)
    ! This subroutine reads precalculated binary obj files (which is faster than reading ascii files)
    implicit none
    character*150 flnm
    character*1 idstr
    integer veccount,facecount,i
    real*4 reado(3)
    integer intdo(3),io,tmp_ctr,tmp_read
    real*4,allocatable :: tmparr(:,:)
    integer, allocatable :: faceid(:)
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old",form="unformatted")
    open(45,file="../Area_IDS.dat",status="old")
    write(*,*) flnm


    read(44) veccount,facecount
    allocate(faceid(facecount))
    tmp_ctr=1
    do

      read(45,*,IOSTAT=io) tmp_read
      faceid(tmp_ctr)=tmp_read
      if(io.lt.0) exit
      tmp_ctr=tmp_ctr+1
    enddo  
    close(45) 
    ctr=tmp_ctr

    allocate(origtria(ctr,10))
    allocate(tmparr(facecount,10))
    read(44) tmparr

    do i=1,ctr
    origtria(i,:)=tmparr(faceid(i),:)
  enddo
    deallocate(tmparr)
    write(*,*) "READ BINOBJFILE COMBE DONE",ctr
  END SUBROUTINE READ_BINOBJFILE_COMBE

  SUBROUTINE READ_BINOBJFILE_ERMAKOV(flnm)
    ! This subroutine reads precalculated binary obj files (which is faster than reading ascii files)
    implicit none
    character*150 flnm
    character*1 idstr
    integer veccount,facecount,i
    real*4 reado(3)
    integer intdo(3),io,tmp_ctr,tmp_read
    real*4,allocatable :: tmparr(:,:)
    integer, allocatable :: faceid(:)
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old",form="unformatted")
    open(45,file="../Area_Ermakov_IDS.dat",status="old")
    write(*,*) flnm


    read(44) veccount,facecount
    allocate(faceid(facecount))
    tmp_ctr=1
    do

      read(45,*,IOSTAT=io) tmp_read
      faceid(tmp_ctr)=tmp_read
      if(io.lt.0) exit
      tmp_ctr=tmp_ctr+1
    enddo  
    close(45) 
    ctr=tmp_ctr

    allocate(origtria(ctr,10))
    allocate(tmparr(facecount,10))
    read(44) tmparr

    do i=1,ctr
    origtria(i,:)=tmparr(faceid(i),:)
  enddo
    deallocate(tmparr)
    write(*,*) "READ BINOBJFILE ERMAKOV DONE",ctr
  END SUBROUTINE READ_BINOBJFILE_ERMAKOV

    SUBROUTINE READ_BINOBJFILE_CUSTOM(flnm,AREA_IDS)
    ! This subroutine reads precalculated binary obj files (which is faster than reading ascii files)
    implicit none
    character*150 flnm,AREA_IDS
    character*1 idstr
    integer veccount,facecount,i
    real*4 reado(3)
    integer intdo(3),io,tmp_ctr,tmp_read
    real*4,allocatable :: tmparr(:,:)
    integer, allocatable :: faceid(:)

    open(44,file=trim(flnm),status="old",form="unformatted")
    open(45,file=AREA_IDS,status="old")
    ! write(*,*) flnm

    read(44) veccount,facecount
    allocate(faceid(facecount))
    tmp_ctr=1
    do

      read(45,*,IOSTAT=io) tmp_read
      faceid(tmp_ctr)=tmp_read
      if(io.lt.0) exit
      tmp_ctr=tmp_ctr+1
    enddo  
    close(45)
    ctr=tmp_ctr

    allocate(origtria(ctr,10))
    allocate(tmparr(facecount,10))
    read(44) tmparr

    do i=1,ctr
    origtria(i,:)=tmparr(faceid(i),:)
  enddo
    deallocate(tmparr)
    write(*,*) "READ BINOBJFILE CUSTOM DONE",ctr
  END SUBROUTINE READ_BINOBJFILE_CUSTOM



    SUBROUTINE READ_BINOBJFILE_LATITUDE(flnm)
    ! This subroutine reads precalculated binary obj files (which is faster than reading ascii files)
    implicit none
    character*150 flnm
    character*1 idstr
    integer veccount,facecount,i
    real*4 reado(3)
    integer intdo(3),io,tmp_ctr,tmp_read
    real*4,allocatable :: tmparr(:,:)
    integer, allocatable :: faceid(:)
    flnm=trim(adjustl(trim(flnm)))
    open(44,file=trim(flnm),status="old",form="unformatted")
    open(45,file="../Area_Latitudinal_IDS.dat",status="old")
    write(*,*) flnm


    read(44) veccount,facecount
    allocate(faceid(facecount))
    tmp_ctr=1
    do

      read(45,*,IOSTAT=io) tmp_read
      faceid(tmp_ctr)=tmp_read
      if(io.lt.0) exit
      tmp_ctr=tmp_ctr+1
    enddo  
    close(45) 
    ctr=tmp_ctr

    allocate(origtria(ctr,10))
    allocate(tmparr(facecount,10))
    read(44) tmparr

    do i=1,ctr
    origtria(i,:)=tmparr(faceid(i),:)
  enddo
    deallocate(tmparr)
    write(*,*) "READ BINOBJFILE LATITUDE DONE",ctr
  END SUBROUTINE READ_BINOBJFILE_LATITUDE

  SUBROUTINE LOAD_METEOROIDS()
    ! This subroutine loads meteoroid files and populates "radmap2"
    implicit none
    real*4 lon,lat,vel,flux,tons(4)
    integer ilon,ilat,ivel,i
    open(21, file="./AST", status="OLD")
    open(22, file="./JFC", status="OLD")
    open(23, file="./HTC", status="OLD")
    open(24, file="./OCC", status="OLD")
    write(*,*) "How many tons of AST,JFC,HTC,OCC (in tons, comma separated)?"
    ! read(*,*) tons
    tons=(/3.7,34.6,2.82,2.18/) ! From Pokorny et al. (2019) - Meteoroids at the Moon in tons/day
    radmap2=0e0
    do i=1,4
      do
        read(20+i,*,end=555) lon,lat,vel,flux 
        ilon=nint((lon+181e0)/2e0)
        ilat=nint((lat+91e0)/2e0)
        ivel=nint((vel+1e0)/2e0)
        radmap2(ilon,ilat,1)=radmap2(ilon,ilat,1)+flux*tons(i)
        radmap2(ilon,ilat,2)=radmap2(ilon,ilat,2)+flux*tons(i)*vel**2*0.5 !ENERGY 1/2*m*v^2
        radmap2(ilon,ilat,3)=radmap2(ilon,ilat,3)+flux*tons(i)*vel**2.46*7.358 ! WITH K&G2001 factor 
      enddo
      555 close(20+i)
    enddo

  END SUBROUTINE LOAD_METEOROIDS

  !   SUBROUTINE LOAD_METEOROIDS_PAR(dir,ifile)
  !   ! This subroutine loads meteoroid files and populates "radmap2"
  !   implicit none
  !   integer,intent(in) :: ifile
  !   real*4 lon,lat,vel,flux,tons(4)
  !   integer ilon,ilat,ivel,i
  !   character*150:: flnm,dir
  !   write(flnm,"(I3)") ifile
  !   flnm=CLEANSPACE(flnm)
  !   dir=CLEANSPACE(dir)
  !   open(21, file=trim(adjustl(trim(dir)))//"AST"//flnm, status="OLD")
  !   open(22, file=trim(adjustl(trim(dir)))//"JFC"//flnm, status="OLD")
  !   open(23, file=trim(adjustl(trim(dir)))//"HTC"//flnm, status="OLD")
  !   open(24, file=trim(adjustl(trim(dir)))//"OCC"//flnm, status="OLD")
  !   ! write(*,*) "How many tons of AST,JFC,HTC,OCC (in tons, comma separated)?"
  !   ! read(*,*) tons
  !   tons=(/3.7,34.6,2.82,2.18/) ! From Pokorny et al. (2019) - Meteoroids at the Moon in tons/day

  !   do i=1,4
  !     do
  !       read(20+i,*,end=555) lon,lat,vel,flux 
  !       ilon=int((lon+181.1e0)/2e0)
  !       ilat=int((lat+91.1e0)/2e0)
  !       ivel=int((vel+1.1e0)/2e0)
  !       radmap2par(ilon,ilat,1,ifile)=radmap2par(ilon,ilat,1,ifile)+flux*tons(i)
  !       radmap2par(ilon,ilat,2,ifile)=radmap2par(ilon,ilat,2,ifile)+flux*tons(i)*vel**2*0.5 !ENERGY 1/2*m*v^2
  !       radmap2par(ilon,ilat,3,ifile)=radmap2par(ilon,ilat,3,ifile)+flux*tons(i)*vel**2.46*7.358 ! WITH K&G2001 factor 
  !     enddo
  !     555 close(20+i)
  !   enddo

  ! END SUBROUTINE LOAD_METEOROIDS_PAR



    SUBROUTINE LOAD_METEOROIDS_PAR_SFD(dir,ifile,sfd)
    ! This subroutine loads meteoroid files and populates "radmap2"
    implicit none
    optional sfd
    integer,intent(in) :: ifile
    real*4 lon,lat,vel,flux,tons(4),area_flux
    integer ilon,ilat,ivel,i
    character*150:: flnm,dir,sfd
    write(flnm,"(I3)") ifile
    flnm=CLEANSPACE(flnm)
    dir=CLEANSPACE(dir)

    IF (PRESENT(sfd)) then
    open(21, file=trim(adjustl(trim(dir)))//"AST"//trim(adjustl(flnm))//"_"//trim(adjustl(sfd)), status="OLD")
    open(22, file=trim(adjustl(trim(dir)))//"JFC"//trim(adjustl(flnm))//"_"//trim(adjustl(sfd)), status="OLD")
    open(23, file=trim(adjustl(trim(dir)))//"HTC"//trim(adjustl(flnm))//"_"//trim(adjustl(sfd)), status="OLD")
    open(24, file=trim(adjustl(trim(dir)))//"OCC"//trim(adjustl(flnm))//"_"//trim(adjustl(sfd)), status="OLD")
    else
    open(21, file=trim(adjustl(trim(dir)))//"AST"//flnm, status="OLD")
    open(22, file=trim(adjustl(trim(dir)))//"JFC"//flnm, status="OLD")
    open(23, file=trim(adjustl(trim(dir)))//"HTC"//flnm, status="OLD")
    open(24, file=trim(adjustl(trim(dir)))//"OCC"//flnm, status="OLD")
    endif  
    ! write(*,*) "How many tons of AST,JFC,HTC,OCC (in tons, comma separated)?"
    ! read(*,*) tons
    tons=(/3.7,34.6,2.82,2.18/) ! From Pokorny et al. (2019) - Meteoroids at the Moon in tons/day

    do i=1,4
      do
        read(20+i,*,end=555) lon,lat,vel,flux,area_flux
        ilon=int((lon+181.1e0)/2e0)
        ilat=int((lat+91.1e0)/2e0)
        ivel=int((vel+1.1e0)/2e0)
        radmap2par(ilon,ilat,1,ifile)=radmap2par(ilon,ilat,1,ifile)+flux*tons(i)
        radmap2par(ilon,ilat,2,ifile)=radmap2par(ilon,ilat,2,ifile)+flux*tons(i)*vel**2*0.5 !ENERGY 1/2*m*v^2
        radmap2par(ilon,ilat,3,ifile)=radmap2par(ilon,ilat,3,ifile)+flux*tons(i)*vel**2.46*7.358 ! WITH K&G2001 factor 
        radmap2par(ilon,ilat,4,ifile)=radmap2par(ilon,ilat,4,ifile)+area_flux*tons(i) ! This is area flux
      enddo
      555 close(20+i)
    enddo

  END SUBROUTINE LOAD_METEOROIDS_PAR_SFD


  SUBROUTINE TIMER()
    implicit none
    ! This little routine tracks time for us
    CALL CPU_TIME(tnow)
    write(*,*) "TOTAL RUNTIME:",tnow,"SEGMENT RUNTIME",tnow-tprev
    tprev=tnow
  END SUBROUTINE TIMER

  SUBROUTINE OPEN_OUTPUT(lonbeg,lonend,latbeg,latend)
    ! Creates an output binary file using the initial parameters
    implicit none
    real*4 lonbeg,lonend,lonstep,latbeg,latend
    character*150 :: flnm2=""

    write(flnm2,"(A6,I4,A1,I4,A1,I4,A1,I4)") "IMPACTS_",nint(lonbeg),"_",nint(lonend),"_",nint(latbeg),"_",nint(latend)
    flnm2=CLEANSPACE(flnm2)
    write(*,*) "The output file is:",flnm2
    open(66,file=trim(flnm2),status="unknown",form="unformatted")
  END SUBROUTINE OPEN_OUTPUT

  SUBROUTINE allocate_and_clear()
    implicit none

  allocate(radmap2par(180,90,4,NFILES))
  radmap2par=0e0

  allocate(radmapAVG(90,4,NFILES))
  radmapAVG=0e0

  allocate(elon(NFILES))
  elon=0d0

  allocate(radmapTOTAL(180,90,4,NFILES))
  radmapTOTAL=0e0

  END SUBROUTINE allocate_and_clear

  SUBROUTINE allocate_object_arrays()
    implicit none
  allocate(ARvhit(comp_ctr,3))
  allocate(ARsurfvec(comp_ctr,3))
  allocate(impactC(comp_ctr))
  allocate(LightMPar(comp_ctr,4,NFILES))

  ARvhit=0.0
  ARsurfvec=0.0
  impactC=0.0
  LightMPar=0.0




end SUBROUTINE allocate_object_arrays

 SUBROUTINE GET_PQ (inc,capom,omega,P,Q)
 real*8 inc,capom,omega
 real*8 P(3),Q(3),DP(3),DQ(3)
  P(1)=dcos(capom)*dcos(omega)-dcos(inc)*dsin(omega)*dsin(capom)
  P(2)=dsin(capom)*dcos(omega)+dcos(inc)*dsin(omega)*dcos(capom)
  P(3)=dsin(inc)*dsin(omega)

  Q(1)=-dcos(capom)*dsin(omega)-dcos(inc)*dcos(omega)*dsin(capom)
  Q(2)=-dsin(capom)*dsin(omega)+dcos(inc)*dcos(omega)*dcos(capom)
  Q(3)=dsin(inc)*dcos(omega)
  end subroutine


SUBROUTINE GET_RVEC(P,Q,a,e,anom,Rvec)
  implicit none
 real*8 P(3),Q(3),a,e,anom,Rvec(3),r
 r=a*(1d0-e*e)/(1d0+e*dcos(anom))
 Rvec=r*(dcos(anom)*P+dsin(anom)*Q)
 end subroutine

 SUBROUTINE GET_VELVEC(P,Q,a_in,e,anom,Vvec)
  implicit none
 real*8 P(3),Q(3),a_in,e,anom,Vvec(3),n,GM
 real*8 a
 a=a_in*1.5e11
 GM=1.32712440018d20
 n=dsqrt(GM*(a)**(-3d0)) !! TO METERS per SECOND
 Vvec=n*a/dsqrt(1d0-e*e)*(-dsin(anom)*P+(e+dcos(anom))*Q)
 end subroutine

 SUBROUTINE LOAD_ORBITAL_ELEMENTS()
  implicit none
  integer i
  real*4 dummy

  open(223,file="orbital_data.in",status="old")
  read(223,*) orbel_vec
  read(223,*) dummy

  allocate(TAA(NFILES))
  TAA=0.0
  do i=1,NFILES
    read(223,*) TAA(i)
  enddo

  TAA=TAA*deg2rad

  do i=3,5
  orbel_vec(i)=orbel_vec(i)*deg2rad
  enddo
  write(*,*) "Orbital elements and true anomalies SUCCESSFULLY loaded!"

END SUBROUTINE


 SUBROUTINE GET_ELONGATIONS()
  implicit none
  real*8 P(3),Q(3)
  real*8 a,e,inc,capom,omega,truean
  real*8 Vvec(3)
  integer i

  write(*,*) "INSIDE_GET_ELONGATIONS"

  write(*,*) TAA
  write(*,*) orbel_vec
  a=orbel_vec(1)
  e=orbel_vec(2)
  !inc=orbel_vec(3)
  inc=0d0 ! This is 0.0 because we are looking for an elongation along the orbit
  capom=orbel_vec(4)
  omega=orbel_vec(5)

  P=0d0
  Q=0d0
  Vvec=0d0
  CALL GET_PQ(inc,capom,omega,P,Q)

  do i=1,NFILES
    truean=TAA(i)
    
    call GET_VELVEC(P,Q,a,e,truean,Vvec)
    elon(i)=atan2(Vvec(2),Vvec(1))*rad2deg-328.249 ! I got this number from ekl2orb.c -it's the ecliptic longitude  along the orbit ! Change this to more general case
    elon(i)=mod(elon(i),360.0)
    if(elon(i).lt.0) elon(i)=elon(i)+360.0
    ! write(*,*) Vvec
    ! write(*,*) i,truean,elon(i)
  END DO

 END SUBROUTINE 

  SUBROUTINE PREPARE_ONE_ORBIT()
    implicit none
    integer i,j,ilon,ilat,i2,jjj,lonii,latii,result(2),lonfini

 do i=1,180
  if(mod(i,18).eq.0) write(*,*) "PREPARING RADIANT MAPS",i/1.8, "% DONE"
    do j=1,90
      ilon=i*2-181
      ilat=j*2-91
  do i2=0,358,2
 

  do jjj=1,NFILES
  result=GET_SPHERE2(ilon*deg2rad,ilat*deg2rad,tilt*deg2rad,0.0,i2) ! Elon = 0.0 because the asteroid does not precess
  
  lonii=result(1)
  latii=result(2)

  lonfini=i+int(elon(jjj)*0.5)
  if(lonfini.lt.1) lonfini=lonfini+180
  if(lonfini.gt.180) lonfini=lonfini-180

  ! write(*,*) lonii,latii,i,j,ilon,ilat,i2,elon(jjj)  
 radmapTOTAL(lonfini,j,1,jjj)=radmapTOTAL(lonfini,j,1,jjj)+radmap2par(lonii,latii,1,jjj)
 radmapTOTAL(lonfini,j,2,jjj)=radmapTOTAL(lonfini,j,2,jjj)+radmap2par(lonii,latii,2,jjj)
 radmapTOTAL(lonfini,j,3,jjj)=radmapTOTAL(lonfini,j,3,jjj)+radmap2par(lonii,latii,3,jjj)
 radmapTOTAL(lonfini,j,4,jjj)=radmapTOTAL(lonfini,j,4,jjj)+radmap2par(lonii,latii,4,jjj)

 enddo

enddo
enddo
enddo

radmapTOTAL=radmapTOTAL/(180.0**1e0) ! This is for i2=0,358,2- simulates the asteroid rotation during one day
  END SUBROUTINE  


  SUBROUTINE PREPARE_AVERAGED_ORBIT()
    implicit none
    integer i,j,ilon,ilat,i2,jjj,lonii,latii,result(2),lonfini,ielon
    real*4 elonloc

do i=1,180
    if(mod(i,18).eq.0) write(*,*) "PREPARING RADIANT MAPS",i/1.8, "% DONE"
    do j=1,90
      ilon=i*2-181
      ilat=j*2-91
  do i2=0,358,2


  do jjj=1,NFILES
  result=GET_SPHERE2(ilon*deg2rad,ilat*deg2rad,tilt*deg2rad,0.0,i2)
  
  lonii=result(1)
  latii=result(2)
  ! write(*,*) lonii,latii,i,j,ilon,ilat,i2,elon(jjj)  
  do ielon=1,359,2
  elonloc=ielon*1.0
  lonfini=i+int(elonloc*0.5)
  if(lonfini.lt.1) lonfini=lonfini+180
  if(lonfini.gt.180) lonfini=lonfini-180

  ! write(*,*) lonii,latii,i,j,ilon,ilat,i2,elon(jjj)  
 radmapTOTAL(lonfini,j,1,jjj)=radmapTOTAL(lonfini,j,1,jjj)+radmap2par(lonii,latii,1,jjj)
 radmapTOTAL(lonfini,j,2,jjj)=radmapTOTAL(lonfini,j,2,jjj)+radmap2par(lonii,latii,2,jjj)
 radmapTOTAL(lonfini,j,3,jjj)=radmapTOTAL(lonfini,j,3,jjj)+radmap2par(lonii,latii,3,jjj)
 radmapTOTAL(lonfini,j,4,jjj)=radmapTOTAL(lonfini,j,4,jjj)+radmap2par(lonii,latii,4,jjj)

 enddo
enddo
enddo
enddo
enddo

radmapTOTAL=radmapTOTAL/(180.0**2e0) ! This is for i2=0,358,2- simulates the asteroid rotation during one day
! This is for ielon=1,359,2 - simulates the asteroid angle change with respect to the meteoroid background 
  END SUBROUTINE

  SUBROUTINE PREPARE_AVERAGED_ORBIT_NOROT()
    implicit none
    integer i,j,ilon,ilat,i2,jjj,lonii,latii,result(2),lonfini,ielon
    real*4 elonloc

do i=1,180
    write(*,*) "PREPARING RADIANT MAPS",i/1.8, "% DONE"
    do j=1,90
      ilon=i*2-181
      ilat=j*2-91

  do jjj=1,NFILES
    ! write(*,*) lonii,latii,i,j,ilon,ilat,i2,elon(jjj)  
 radmapTOTAL(i,j,1,jjj)=radmapTOTAL(i,j,1,jjj)+radmap2par(i,j,1,jjj)
 radmapTOTAL(i,j,2,jjj)=radmapTOTAL(i,j,2,jjj)+radmap2par(i,j,2,jjj)
 radmapTOTAL(i,j,3,jjj)=radmapTOTAL(i,j,3,jjj)+radmap2par(i,j,3,jjj)
 radmapTOTAL(i,j,4,jjj)=radmapTOTAL(i,j,4,jjj)+radmap2par(i,j,4,jjj)

enddo
enddo
enddo

! radmapTOTAL=radmapTOTAL/(180.0**2e0) ! This is for i2=0,358,2- simulates the asteroid rotation during one day
! This is for ielon=1,359,2 - simulates the asteroid angle change with respect to the meteoroid background 
  END SUBROUTINE 

    SUBROUTINE PREPARE_AVERAGED_ORBIT_MERCURY_2ORB()
    implicit none
    integer i,j,ilon,ilat,i2,jjj,lonii,latii,result(2),lonfini,ielon
    real*4 elonloc

do i=1,180
    write(*,*) "PREPARING RADIANT MAPS",i/1.8, "% DONE"
    do j=1,90
      ilon=i*2-181
      ilat=j*2-91

  do jjj=1,NFILES
    ! write(*,*) lonii,latii,i,j,ilon,ilat,i2,elon(jjj)  
 radmapTOTAL(i,j,1,jjj)=radmapTOTAL(i,j,1,jjj)+radmap2par(i,j,1,jjj)
 radmapTOTAL(i,j,2,jjj)=radmapTOTAL(i,j,2,jjj)+radmap2par(i,j,2,jjj)
 radmapTOTAL(i,j,3,jjj)=radmapTOTAL(i,j,3,jjj)+radmap2par(i,j,3,jjj)
 radmapTOTAL(i,j,4,jjj)=radmapTOTAL(i,j,4,jjj)+radmap2par(i,j,4,jjj)

 ilon=181-i
 radmapTOTAL(ilon,j,1,jjj)=radmapTOTAL(ilon,j,1,jjj)+radmap2par(i,j,1,jjj)
 radmapTOTAL(ilon,j,2,jjj)=radmapTOTAL(ilon,j,2,jjj)+radmap2par(i,j,2,jjj)
 radmapTOTAL(ilon,j,3,jjj)=radmapTOTAL(ilon,j,3,jjj)+radmap2par(i,j,3,jjj)
 radmapTOTAL(ilon,j,4,jjj)=radmapTOTAL(ilon,j,4,jjj)+radmap2par(i,j,4,jjj)
enddo
enddo
enddo

radmapTOTAL=radmapTOTAL*0.5d0 !Average over 2 orbits
! This is for ielon=1,359,2 - simulates the asteroid angle change with respect to the meteoroid background 
  END SUBROUTINE   


 SUBROUTINE ERROR_INPUT_1()
  implicit none
  write(*,*) "Your input is different from permitted inputs 1, 2, and 3. Program will be terminated"
  write(*,*)
  stop
END SUBROUTINE ERROR_INPUT_1

 SUBROUTINE ERROR_INPUT_2()
  implicit none
  write(*,*) "Your input is different from permitted inputs 1, and 2 Program will be terminated"
  write(*,*)
  stop
END SUBROUTINE ERROR_INPUT_2




end module my_functions
 
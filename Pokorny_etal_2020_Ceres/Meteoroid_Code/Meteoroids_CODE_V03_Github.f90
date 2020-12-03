! Meteoroids ray-tracing code for OBJ files using FASTBVH and TINYOBJ
! Written by Petr Pokorny (petr.pokorny@nasa.gov / pokorny@cua.edu) 
!
! Disclaimer: This code is not final (but works) and will be properly commented soon.
! V03 - Dec 1st, 2020

program test
  use my_functions
  use omp_lib
  implicit none
  character(len=*), parameter :: obj_key = "CERES"
  character(len=*), parameter :: path_obj = "../../../ICQtoOBJ/CERES_0512.obj"
  character(len=150), parameter :: path_bin = "../../../ICQtoOBJ/CERES_0512.bin"
  character*150 :: flnm="",dir="", sfd="",AREA_IDS
  character*3 :: chnm
  ! real, dimension(3) :: rayorig = [0, 0, 10000]
  real, dimension(3) :: vhit = [0, 0, -1],v1,v2,v3,rayvec,v1v2,v1v3
  real surfvec(3),centroid(3)
  real :: cosdec,prec
  integer i,i1,j,ilon,ilat,lonii,latii,jjj,i2,result(2),ielon
  integer :: strlon,endlon,strlat,endlat,stp,lonfini,iprec,selector
  integer :: select_orbit, select_area


  call omp_set_num_threads(4) ! Sets the number of threads used in the code - this code is memory heavy so 4 threads worked the best for me
  write(*,*)'Using the following number of threads for this code:',omp_get_max_threads(); write(*,*)





  write(*,*) "What is the starting longitude, final longitude, starting latitude,&
              ending latitude, step in lat/lon, axial tilt in degrees?"
  read(*,*) strlon,endlon,strlat,endlat,stp,tilt

  write(*,*) "Write '1' if you want the one orbit calculation"
  write(*,*) "Write '2' if you want the averaged orbit calculation"
  write(*,*) "Write '3' if you want the 3:2 spin-orbit calculation"
  read(*,*) select_orbit
  if(((select_orbit.eq.1).or.(select_orbit.eq.2).or.(select_orbit.eq.3)).eqv..false.) CALL ERROR_INPUT_1()


  if(select_orbit.eq.1) dir="../METEOROIDS_DATA_169/"
  if(select_orbit.eq.2) dir="../METEOROIDS_DATA_AVG/"

  write(*,*) "Write '1' if you want to analyze the entire asteroid"
  write(*,*) "Write '2' if you want to analyze only selected faces"

  read(*,*) select_area
  if(((select_area.eq.1).or.(select_area.eq.2)).eqv..false.) CALL ERROR_INPUT_2()

  write(*,*) "Number of analyzed true anomaly angles (for averaged orbits write 1)"
  read(*,*) NFILES


  call allocate_and_clear()
  call LOAD_ORBITAL_ELEMENTS()
  call get_elongations()

  write(*,*) "Size-frequency distribution index (write e.g., 4.0) - this is a part of the name of the meteoroid input file"
  read(*,*)SFD

  do jjj=1,NFILES
    call LOAD_METEOROIDS_PAR_SFD(dir,jjj,sfd) !!! FILLS RAD2MAP with numbers
  write(*,*) "Loaded Meteoroid Input Files Number:",jjj,"with SFD:",sfd
  enddo


   write(*,*) "------------------------------"
   write(*,*) "Settings:",strlon,endlon,strlat,endlat,stp,tilt
   write(*,*) "Selectors: Orbit: ",select_orbit,", Area: ",select_area
   write(*,*) "Number of processed true anomaly files:",NFILES
   



  if(select_area.eq.1) CALL READ_BINOBJFILE(path_bin)
  if(select_area.eq.2) then 
    read(*,*) AREA_IDS
  write(*,*) "IDs of the analyzed areas:",AREA_IDS
    call READ_BINOBJFILE_CUSTOM(path_bin,AREA_IDS)
  endif

   write(*,*) "Size-frequency distribution index:",SFD
   write(*,*) "------------------------------"
   write(*,*) 


   if(select_orbit.eq.1) CALL PREPARE_ONE_ORBIT()
   if(select_orbit.eq.2) CALL PREPARE_AVERAGED_ORBIT()
   if(select_orbit.eq.3) CALL PREPARE_AVERAGED_ORBIT_MERCURY_2ORB()


  call load_obj(obj_key, path_obj)
  write(*,"(A)",advance = 'no') "EVERYTHING LOADED - "; call timer()

  
  
do jjj=1,NFILES
  write(chnm,"(I3)") jjj
  chnm=adjustl(trim(adjustl(chnm)))
  open(2222,file="RADIANTS_"//trim(chnm)//".txt",status='unknown')
  open(2223,file="RADIANTS_PAR_"//trim(chnm)//".txt",status='unknown')
    do i=1,180
    do j=1,90
  write(2222,*) i,j,radmapTOTAL(i,j,1,jjj),radmapTOTAL(i,j,2,jjj),radmapTOTAL(i,j,3,jjj),radmapTOTAL(i,j,4,jjj)
  write(2223,*) i,j,radmap2par(i,j,1,jjj),radmap2par(i,j,2,jjj),radmap2par(i,j,3,jjj),radmap2par(i,j,4,jjj)
 enddo
enddo
enddo
!!! TEST 

  comp_ctr = ctr
  call allocate_object_arrays()


    ! call timer()
        do i1=1,ctr
        v1=origtria(i1,1:3) ! We are going through all triangles based on their rotated z-axis distance
        v2=origtria(i1,4:6)
        v3=origtria(i1,7:9)
        vhit=real(POINT_INSIDE_TRIANGLE((v1),(v2),(v3)))
        ARvhit(i1,1:3)=vhit
        v1v2=v2-v1
        v1v3=v3-v1
        surfvec=cross(v1v2,v1v3)
        if(SCALAR(surfvec,v1).lt.0) surfvec=-surfvec
        surfvec=surfvec/sqrt(scalar(surfvec,surfvec))
        ARsurfvec(i1,1:3)=surfvec
      enddo

      do ilon=strlon,endlon,stp
        do ilat=strlat,endlat,stp

          impactC=0e0
          rayvec=GET_RAYVEC(real(ilon)*deg2rad,ilat*deg2rad)

          !$omp parallel &
          !$omp shared (impactC,ctr,ARvhit,rayvec,ARsurfvec) &       
          !$omp private (i1,cosdec)

          !$omp do reduction ( + : impactC )
          do i1=1,ctr
            cosdec=SCALAR(rayvec,ARsurfvec(i1,:))
            if(cosdec.gt.0d0) then
              if((query_occluded(obj_key, ARvhit(i1,:)+rayvec*0.01, rayvec)).eqv.(.true.))  THEN

              endif
            else
              cosdec=0e0
            endif
            impactC(i1)=cosdec  
          enddo
          !$omp end do

          !$omp end parallel

            write(*,"(A,I4,1X,A,I4)",advance='no') "Longitude: ",ilon,", Latitude: ", ilat; call timer()

          lonii=nint((ilon+181)*0.5)
          latii=nint((ilat+91)*0.5)
          !$omp parallel &
          !$omp shared (impactC,ctr,radmapTOTAL,LightMPar,tilt,lonii,latii) &       
          !$omp private (i,jjj)
          !$omp do reduction ( + : LightMPar )

          do jjj=1,NFILES
           do i=1,ctr
            LightMPar(i,1,jjj)=LightMPar(i,1,jjj)+impactC(i)*radmapTOTAL(lonii,latii,1,jjj)
            LightMPar(i,2,jjj)=LightMPar(i,2,jjj)+impactC(i)*radmapTOTAL(lonii,latii,2,jjj)
            LightMPar(i,3,jjj)=LightMPar(i,3,jjj)+impactC(i)**3e0*radmapTOTAL(lonii,latii,3,jjj)
            LightMPar(i,4,jjj)=LightMPar(i,4,jjj)+impactC(i)*radmapTOTAL(lonii,latii,4,jjj)
          enddo
        enddo

        !$omp end do

        !$omp end parallel


      enddo
    enddo


do jjj=1,NFILES
  write(flnm,"(I3)") jjj
  flnm=CLEANSPACE(flnm)
  open(221,file="METEOROIDS_BIN_"//flnm,status="unknown",form="unformatted")
  do i=1,ctr
    centroid = POINT_INSIDE_TRIANGLE(origtria(i,1:3),origtria(i,4:6),origtria(i,7:9))
    write(221) real(centroid),real(LightMPar(i,1:4,jjj))*(stp*0.5e0)**2e0
  enddo

enddo
  call free_obj(obj_key)

  contains

  ! ====================================
! C BINDING FUNCTIONS AND SUBROUTINES 
! ====================================

  subroutine load_obj(obj_key, path)
    use iso_c_binding, only: c_null_char
    
    character(len=*) :: obj_key
    character(len=*) :: path

    interface
       subroutine sv_load_obj(obj_key, path)
         use iso_c_binding, only: c_char
         character(kind=c_char, len=*) :: obj_key
         character(kind=c_char, len=*) :: path
       end subroutine sv_load_obj
    end interface

    call sv_load_obj(obj_key // c_null_char, path // c_null_char)
  end subroutine load_obj

  subroutine free_obj(obj_key)
    use iso_c_binding, only: c_char, c_null_char

    character(len=*) :: obj_key

    interface
       subroutine sv_free_obj(obj_key)
         use iso_c_binding, only: c_char
         character(kind=c_char, len=*) :: obj_key
       end subroutine sv_free_obj
    end interface

    call sv_free_obj(obj_key // c_null_char)
  end subroutine free_obj

  logical function query_occluded(obj_key, r, d)
    use iso_c_binding, only: c_bool, c_char, c_null_char

    character(len=*) :: obj_key
    real, dimension(3) :: r
    real, dimension(3) :: d

    interface
       logical(c_bool) function sv_query_occluded(obj_key, r, d)
         use iso_c_binding, only: c_bool, c_char
         character(kind=c_char, len=*) :: obj_key
         real, dimension(*) :: r
         real, dimension(*) :: d
       end function sv_query_occluded
    end interface

    query_occluded = sv_query_occluded(obj_key // c_null_char, r, d)
  end function query_occluded
  
end program test

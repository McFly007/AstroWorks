 PROGRAM OBJ2BIN
 	implicit none
 	character*1 idstr
 	real*4 reado(3),vecarr(10000000,3),facearr(10000000,3),origtria(10000000,10)
 	integer :: intdo(3),facecount,veccount
 	open(44,file="CERES_0512.obj",status="old")
	do 
    read(44,*,end=666) idstr,reado(1),reado(2),reado(3)
    ! write(*,*) idstr,reado(1),reado(2),reado(3)
    if(idstr.eq."v") then
      veccount=veccount+1
      vecarr(veccount,:)=reado
    elseif(idstr.eq."f") then
      intdo(1)=nint(reado(1))
      intdo(2)=nint(reado(2))
      intdo(3)=nint(reado(3))
      facecount=facecount+1
      ! facearr(facecount,:)=intdo
      origtria(facecount,:)=(/vecarr(intdo(1),:),vecarr(intdo(2),:),vecarr(intdo(3),:),0e0/)
    endif  
 	enddo
 	666 close(44)
 	open(45,file="CERES_0512.bin",status="unknown",form="unformatted")
 	write(45) veccount,facecount
 	! write(45) vecarr(1:veccount,:)
 	! write(45) facearr(1:facecount,:)
 	write(45) origtria(1:facecount,:)
 	close(45)

 END


 MODULE my_module
 	implicit none
	integer, allocatable :: start_p(:),end_p(:),pivots(:)
  integer :: NAREA

	contains
	FUNCTION GIVE_PIVOT(i)
	implicit none
	integer i,GIVE_PIVOT,kk
	
	do kk=1,NAREA
	if((i.ge.start_p(kk)).and.(i.le.end_p(kk))) GIVE_PIVOT=kk
	enddo
	! write(*,*) i,GIVE_PIVOT
END FUNCTION GIVE_PIVOT



  recursive subroutine quicksort(a,col)
! This is a Quicksort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.)
! I modified it to arrays with multiple columns
    implicit none
    real*4 :: a(:,:) ! Array we want to sort
    integer :: col ! Column that will control the sorting
    real*4 x, t(size(a,dim=2))
    integer :: first = 1, last
    integer i, j

    ! write(*,*) size(a,dim=2),size(t)
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


	END MODULE

 PROGRAM ANALYZE_COMBE
 	use my_module
 	implicit none

 	real*4, allocatable ::  x(:),y(:),z(:),f1(:),f2(:),f3(:),f4(:),tmp_arr(:,:)
 	real*4, allocatable :: sumx(:,:),sumxx(:,:),mean(:,:),variance(:,:),stddev(:,:)
 	real*4 resvec(20)
 	character*10 fl_num,out_num
  character*150 PIVOT_FILE
 	integer NFILES,i,io,n,j,k,i1,i2,l,res_id,pivot

  write(*,*) "Number of files to process: NFILES"
  read(*,*) NFILES


  write(*,*) "Pivot file to process"
  read(*,*) PIVOT_FILE
  open(21,file=trim(adjustl(PIVOT_FILE)),status="old")
  i=1
  do
  read(21,*,end=345) pivot
  i=i+1
  enddo

345 NAREA=i-1

close(21)
allocate(pivots(NAREA))
allocate(start_p(NAREA))
allocate(end_p(NAREA))
allocate(sumx(7,NAREA))
allocate(sumxx(7,NAREA))
allocate(mean(7,NAREA))
allocate(variance(7,NAREA))
allocate(stddev(7,NAREA))

 	open(21,file=trim(adjustl(PIVOT_FILE)),status="old")
 	do i=1,NAREA
 	read(21,*,end=345) pivots(i)
 	if(i.eq.1) then
 		start_p(i)=1
 		end_p(i)=pivots(i)
 	endif
 	if(i.gt.1) then
 	start_p(i)=start_p(i-1)+pivots(i-1)
 	end_p(i)=end_p(i-1)+pivots(i)
 	endif
 	enddo

write(*,*) "PROCESSING -  NAREA:",NAREA,", NFILES:",NFILES

        do j=1,NAREA
          write(out_num,"(I10)") j
          open(321,file="Area_Results_"//trim(adjustl(out_num))//".txt",status="unknown")
          write(321,*) "# Area",j," from Combe et al. (2019) results"
          close(321)
        enddo
 	allocate(x(end_p(NAREA)))
 	allocate(y(end_p(NAREA)))
 	allocate(z(end_p(NAREA)))
 	allocate(f1(end_p(NAREA)))
 	allocate(f2(end_p(NAREA)))
 	allocate(f3(end_p(NAREA)))
  allocate(f4(end_p(NAREA)))



 	do i=1,NFILES
 		write(fl_num,"(I10)") i
 		write(*,*) "Processing file:",i
 		open(22,file="METEOROIDS_BIN_"//trim(adjustl(fl_num)), status="OLD",form="unformatted")

 		x=0.0
 		y=0.0
 		z=0.0
 		f1=0.0
 		f2=0.0
 		f3=0.0

 		do j=1,end_p(NAREA)
 		read(22,IOSTAT=io) x(j),y(j),z(j),f1(j),f2(j),f3(j),f4(j)
 		! write(*,*) "METEOROIDS_BIN_"//trim(adjustl(fl_num)), x(j),y(j),z(j),f1(j),f2(j),f3(j),j,io
    ! read(*,*)
 		if(io.lt.0) exit
 		enddo
 		
 		sumx=0.0
 		sumxx=0.0
 		! NOW I READ ONE FILE - DO SOME STATS
 		do j=1,end_p(NAREA)
 			k=GIVE_PIVOT(j)
 			sumx(:,k)=sumx(:,k)+(/x(j),y(j),z(j),f1(j),f2(j),f3(j),f4(j)/)
 			sumxx(:,k)=sumxx(:,k)+(/x(j)**2.0,y(j)**2.0,z(j)**2.0,f1(j)**2.0,f2(j)**2.0,f3(j)**2.0,f4(j)**2.0/)
 		enddo

        do j=1,NAREA
          write(*,*) "PROCESSING AREA:",j
        	write(out_num,"(I10)") j
        	open(321,file="Area_Results_"//trim(adjustl(out_num))//".txt",position="append")
        n=end_p(j)-start_p(j)+1
 	    	Mean = sumx / (n*1.0)
        Variance = (sumxx - sumx*sumx/(n*1.0))/(n*1.0-1.0)
        StdDev   = SQRT(Variance)	
        ! write(*,*) "Area",j,mean(1,j),stddev(1,j),minval(x(start_p(j):end_p(j))),maxval(x(start_p(j):end_p(j))),&
        ! start_p(j),end_p(j)

        i1=start_p(j)
        i2=end_p(j)
        ! write(*,*) i1,i2
        allocate(tmp_arr((i2-i1+1),7))

        do l=i1,i2
        tmp_arr(l+1-i1,:)=(/x(l),y(l),z(l),f1(l),f2(l),f3(l),f4(l)/)
		enddo

   		res_id=4
   		! write(*,*) size(tmp_arr)
   		call quicksort(tmp_arr,res_id)

   		resvec(1:5)=(/mean(res_id,j),stddev(res_id,j),&
   			tmp_arr(nint((i2-i1+1)*0.05),res_id),&
   			tmp_arr(nint((i2-i1+1)*0.50),res_id),&
 			tmp_arr(nint((i2-i1+1)*0.95),res_id)/)

   		   		res_id=5
   		call quicksort(tmp_arr,res_id)

   		resvec(6:10)=(/mean(res_id,j),stddev(res_id,j),&
   			tmp_arr(nint((i2-i1+1)*0.05),res_id),&
   			tmp_arr(nint((i2-i1+1)*0.50),res_id),&
 			tmp_arr(nint((i2-i1+1)*0.95),res_id)/)

   		   		res_id=6
   		call quicksort(tmp_arr,res_id)

   		resvec(11:15)=(/mean(res_id,j),stddev(res_id,j),&
   			tmp_arr(nint((i2-i1+1)*0.05),res_id),&
   			tmp_arr(nint((i2-i1+1)*0.50),res_id),&
 			tmp_arr(nint((i2-i1+1)*0.95),res_id)/)

            res_id=7
      call quicksort(tmp_arr,res_id)

      resvec(16:20)=(/mean(res_id,j),stddev(res_id,j),&
        tmp_arr(nint((i2-i1+1)*0.05),res_id),&
        tmp_arr(nint((i2-i1+1)*0.50),res_id),&
      tmp_arr(nint((i2-i1+1)*0.95),res_id)/)

        deallocate(tmp_arr)


        write(321,*) resvec
        close(321)
 
    enddo
 		close(22)





 	enddo

 end

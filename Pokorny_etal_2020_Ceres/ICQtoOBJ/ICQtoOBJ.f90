  PROGRAM ICQ2OBJ
    implicit none
    integer q,f,j,i,n0,k
    real*8,allocatable :: v(:,:,:,:)
    integer, allocatable :: n(:,:,:)
    real*8 :: w1(3),w2(3),w3(3),z1,z2

    open(10,file="CERES_SPC181019_0512.ICQ",status="OLD")
    open(20,file="CERES_0512.obj")
    open(21,file="CERES_0512.info")

       read(10,*) q

       allocate(v(3,q+1,q+1,6))
       allocate(n(q+1,q+1,6))
        write(21,*) 6*(q+1)**2, 12*q**2   ! number of vertices and facets
        n0=0      
        do f=1,6
        do j = 0+1,q+1
        do i = 0+1,q+1
          ! write(*,*) k,i,j,f
          read(10,*) (v(k,i,j,f), k=1,3)
          n0=n0+1
          write(20,fmt='(A1,1x,3(f11.5))') "v", (v(k,i,j,f), k=1,3)
          n(i,j,f)=n0
        enddo
        enddo
        enddo
        do i=1+1,q-1+1
          n(i+1,q+1,6)=n(q-i+1,q+1,4)
          n(i+1,0+1,6)=n(i+1,q+1,2)
          n(i+1,0+1,5)=n(q+1,q-i+1,1)
          n(i+1,0+1,4)=n(q-i+1,0+1,1)
          n(i+1,0+1,3)=n(0+1,i+1,1)
          n(i+1,0+1,2)=n(i+1,q+1,1)
        enddo
        do j=1,q-1
          n(q+1,j+1,6)=n(j+1,q+1,5)
          n(q+1,j+1,5)=n(0+1,j+1,4)
          n(q+1,j+1,4)=n(0+1,j+1,3)
          n(q+1,j+1,3)=n(0+1,j+1,2)
          n(0+1,j+1,6)=n(q-j+1,q+1,3)
          n(0+1,j+1,5)=n(q+1,j+1,2)
        enddo
        n(0+1,0+1,3)=n(0+1,0+1,1)
        n(q+1,0+1,4)=n(0+1,0+1,1)
        n(0+1,0+1,2)=n(0+1,q+1,1)
        n(q+1,0+1,3)=n(0+1,q+1,1)
        n(0+1,0+1,4)=n(q+1,0+1,1)
        n(q+1,0+1,5)=n(q+1,0+1,1)
        n(0+1,0+1,5)=n(q+1,q+1,1)
        n(q+1,0+1,2)=n(q+1,q+1,1)
        n(0+1,0+1,6)=n(0+1,q+1,2)
        n(q+1,q+1,3)=n(0+1,q+1,2)
        n(0+1,q+1,5)=n(q+1,q+1,2)
        n(q+1,0+1,6)=n(q+1,q+1,2)
        n(q+1,q+1,4)=n(0+1,q+1,3)
        n(0+1,q+1,6)=n(0+1,q+1,3)
        n(q+1,q+1,5)=n(0+1,q+1,4)
        n(q+1,q+1,6)=n(0+1,q+1,4)
        n0=0
        do f=1,6
        do i=0+1,q-1+1
        do j=0+1,q-1+1
          w1(1)=v(2,i,j,f)*v(3,i+1,j+1,f)&
           -v(3,i,j,f)*v(2,i+1,j+1,f)
          w1(2)=v(3,i,j,f)*v(1,i+1,j+1,f)&
           -v(1,i,j,f)*v(3,i+1,j+1,f)
          w1(3)=v(1,i,j,f)*v(2,i+1,j+1,f)&
           -v(2,i,j,f)*v(1,i+1,j+1,f)
          w2(1)=v(2,i+1,j,f)*v(3,i,j+1,f)&
           -v(3,i+1,j,f)*v(2,i,j+1,f)
          w2(2)=v(3,i+1,j,f)*v(1,i,j+1,f)&
           -v(1,i+1,j,f)*v(3,i,j+1,f)
          w2(3)=v(1,i+1,j,f)*v(2,i,j+1,f)&
           -v(2,i+1,j,f)*v(1,i,j+1,f)
          z1=w1(1)**2+w1(2)**2+w1(3)**2
          z2=w2(1)**2+w2(2)**2+w2(3)**2
          if(z1.le.z2) then
            n0=n0+1
            write(20,fmt=1000) "f", n(i,j,f), &
                               n(i+1,j+1,f), n(i+1,j,f)
            n0=n0+1
            write(20,fmt=1000) "f", n(i,j,f), &
                               n(i,j+1,f), n(i+1,j+1,f)
          else
            n0=n0+1
            write(20,fmt=1000) "f", n(i,j,f), &
                               n(i,j+1,f), n(i+1,j,f)
            n0=n0+1
            write(20,fmt=1000) "f", n(i+1,j,f), &
                               n(i,j+1,f), n(i+1,j+1,f)
          endif
        enddo
        enddo
        enddo
1000 FORMAT(A1,1X,3I10)
      END
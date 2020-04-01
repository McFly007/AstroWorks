   integer(16) :: r_num
   parameter(r_num=400)
   real*8 r_boxminX
   real*8 r_boxmaxX
   real*8 r_boxminY
   real*8 r_boxmaxY
   real*8 :: r_dx
   real*8 :: r_dy
   integer(16) :: r_pivots(r_num*r_num)=0
   integer(16) :: r_numbs(r_num*r_num)=0
   real*8,allocatable :: r_tria(:,:) 
   integer(16) :: r_id(4)=0
   integer(16) :: l1,l2,r_index

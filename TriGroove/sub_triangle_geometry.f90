subroutine triangle_geometry
	use constant
	use global
	implicit none

!	integer lx, ly													! 阶梯的宽度和高度
	integer num_ladder												! 近似三角边的阶梯数
	integer num_at													! 子计算区域总数

	integer i, j, k, ipre, jpre, inum_a
	integer inum_a1, inum_a2, inum_a3

!	if (h_groove<=w_groove) then		! h_groove与w_groove必须为整数倍关系
!		ly = 1
!		lx = ly *(w_groove/h_groove)
!	else
!		lx = 1
!		ly = lx *(h_groove/w_groove)
!	end if

	num_ladder = h_groove/ly
	num_at = 3*num_ladder +1

	if (num_at /= num_a) then
		write(*,*) "=================>error!"
		stop
	end if

	do i=1,num_ladder
		inum_a1 = (i-1)*3 +1
		inum_a2 = (i-1)*3 +2
		inum_a3 = (i-1)*3 +3
		
		if (i==1) then
			gcb(inum_a1, 1) = 0
			gcb(inum_a1, 2) = w_wall
			gcb(inum_a1, 3) = h_groove
			gcb(inum_a1, 4) = h_groove +h_plasma

			IB(inum_a1, 1) = 4
			IB(inum_a1, 2) = 3
			IB(inum_a1, 3) = 2
			IB(inum_a1, 4) = 5

			gcb(inum_a2, 1) = w_wall
			gcb(inum_a2, 2) = w_wall +2*w_groove
			gcb(inum_a2, 3) = h_groove
			gcb(inum_a2, 4) = h_groove +h_plasma

			IB(inum_a2, 1) = 3
			IB(inum_a2, 2) = 3
			IB(inum_a2, 3) = 3
			IB(inum_a2, 4) = 5

			gcb(inum_a3, 1) = w_wall +2*w_groove
			gcb(inum_a3, 2) = w_wall +2*w_groove +w_wall
			gcb(inum_a3, 3) = h_groove
			gcb(inum_a3, 4) = h_groove +h_plasma

			IB(inum_a3, 1) = 3
			IB(inum_a3, 2) = 4
			IB(inum_a3, 3) = 2
			IB(inum_a3, 4) = 5
		else
			gcb(inum_a1, 1) = w_wall +lx*(i-2)
			gcb(inum_a1, 2) = w_wall +lx*(i-1)
			gcb(inum_a1, 3) = h_groove -ly*(i-1)
			gcb(inum_a1, 4) = h_groove -ly*(i-2)

			IB(inum_a1, 1) = 2
			IB(inum_a1, 2) = 3
			IB(inum_a1, 3) = 2
			IB(inum_a1, 4) = 3

			gcb(inum_a2, 1) = w_wall +lx*(i-1)
			gcb(inum_a2, 2) = w_wall +2*w_groove -lx*(i-1)
			gcb(inum_a2, 3) = h_groove -ly*(i-1)
			gcb(inum_a2, 4) = h_groove -ly*(i-2)

			IB(inum_a2, 1) = 3
			IB(inum_a2, 2) = 3
			IB(inum_a2, 3) = 3
			IB(inum_a2, 4) = 3

			gcb(inum_a3, 1) = w_wall +2*w_groove -lx*(i-1)
			gcb(inum_a3, 2) = w_wall +2*w_groove -lx*(i-2)
			gcb(inum_a3, 3) = h_groove -ly*(i-1)
			gcb(inum_a3, 4) = h_groove -ly*(i-2)

			IB(inum_a3, 1) = 3
			IB(inum_a3, 2) = 2
			IB(inum_a3, 3) = 2
			IB(inum_a3, 4) = 3
		end if

	end do

	gcb(num_a, 1) = w_wall +lx*(num_ladder+1-2)
	gcb(num_a, 2) = w_wall +2*w_groove -lx*(num_ladder+1-2)
	gcb(num_a, 3) = h_groove -ly*(num_ladder+1-1)
	gcb(num_a, 4) = h_groove -ly*(num_ladder+1-2)

	IB(num_a, 1) = 2
	IB(num_a, 2) = 2
	IB(num_a, 3) = 2
	IB(num_a, 4) = 3

	gcb = gcb * cells_per_lmdd
	CB = gcb * lmdd0/2.

	NZ(:) = gcb(:, 2) -gcb(:, 1)
	NR(:) = gcb(:, 4) -gcb(:, 3)

	z1	=	0.
	z2	=	(w_groove+w_wall)*lmdd0*2.
	r0	=	0.021-h_groove*lmdd0
	r1	=	0.
	r2	=	(h_groove+h_plasma)*lmdd0

	nzt	=	2*(w_groove+w_wall) * cells_per_lmdd
	nrt	=	(h_groove +h_plasma) * cells_per_lmdd

	if (nzt/=imax1 .OR. nrt/=imax2) then
		write(*,*) "===========================>error!"
		stop
	end if


	if (NZ(num_a) /= NZ(num_a-2)) then
		write(*,*) "======================> Error"
		stop
	end if

	write(*,*) "======================> Success"

end subroutine triangle_geometry
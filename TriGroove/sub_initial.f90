! ++++++++++++++++++++++++++++++++++++++++++++++
!
! 2009-05-08：对获得cmprg_gradient的代码进行了改善，更正了一个错误；
! 
! 2009-05-15 更新：
! 求解Poisson方程时，根据节点位置的不同要使用不同的离散差分格式（特别对于边界节点），因此需要对节点类型和位置做出判断。
! 由如下三个变量的组合来进行判断：
!
! (1) ro_boundary: 标示节点类型.
! =1: 内节点；=2: 边界节点（非角点）；=4: 内角点；=4/3: 外角点
!
! (2) cmprg_gradient：标示边界节点（非角点）处于哪个边界（左、右、上、下），以及角点是由哪两条边界构成.
! =2，左边界；=4，右边界；=5，下边界；=6，上边界。
! 注：1) 由于角点是两个边界的公共点，因此cmprg_gradient的数值要叠加；
!     2) interior边界上的节点相当于内点，其cmprg_gradient均等于1。
!
! (3) cmprg_bdtype：标示边界类型
! cmprg_bdtype值与ib值相同，在角点处亦要叠加，具体叠加方法见代码相关部分；
! 对于内角点和interior边界，cmprg_bdtype = 3。
!
! ++++++++++++++++++++++++++++++++++++++++++++++

SUBROUTINE INITIAL
	use constant
	use global
	use ion
	use grid
	use grid_para
	use magnetic
	use diag
	use mcc_para
	use neutral
	use discharge
	use mod_dadi
	use wall
	use poisson
	
	implicit none
	
	real :: RANUM = 0.95
	real ra,rb

	real lcb(num_a, 4), volumet(num_a), bt(0:IMAX1,0:IMAX2) ! lcb是每条边的长度
	real fluxrho_i, flux
	real krypton(energy_max,4) ! 氪气的碰撞截面，第2，3，4列分别代表弹性、激发和电离碰撞截面
	real, allocatable :: B(:,:)
	integer i, j, ii, jj, inum_a, k, insp, ipre, jpre, num_ea(num_a)
	integer inject_num, fnum_inject

	dzr = (z2-z1) / nzt ! 网格的横向尺寸
	drr = (r2-r1) / nrt ! 网格的径向尺寸
	do i = 1,num_a
		dz(i) = (cb(i,2)-cb(i,1)) / nz(i) ! 第i个子计算区域的网格的横向尺寸
		dr(i) = (cb(i,4)-cb(i,3)) / nr(i) ! 第i个子计算区域的网格的径向尺寸
	end do
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!                   lcb(i,4)
	!             _________________
	!            |                 |
	!  lcb(i,1)  |    REGION i     | lcb(i,2)
	!            |_________________|
	!
	!                  lcb(i,3)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i = 1,num_a
		lcb(i,1) = cb(i,4) - cb(i,3)
		lcb(i,2) = cb(i,4) - cb(i,3)
		lcb(i,3) = cb(i,2) - cb(i,1)
		lcb(i,4) = cb(i,2) - cb(i,1)
	end do
	
	! ro_boundary：计算除了interior类型之外的边界处的密度时需要的系数，同时作为计算Possion方程的判断条件之一
	! 一般地：ro_boundary = 2;
	! 特别地，对于边界交点，ro_boundary = 4 or 4/3.
	ro_boundary = 0.
	do inum_a = 1,num_a
	    ! 获得每个子计算区域的四条边上节点（不包括边界交点）处的ro_boundary
	    do i = 1,2
		    jpre = gcb(inum_a,i+2)      ! jpre：计算区域inum_a的上、下边界的行索引
		    if(ib(inum_a,i+2) /= 3) then        ! 非interior边界
		        do ipre = gcb(inum_a,1)+1,gcb(inum_a,2)-1       ! ipre：jpre边界上的节点（不包括边界端点）
		            ro_boundary(ipre,jpre) = 2.
		        end do
		    end if
		    ipre = gcb(inum_a,i)        ! ipre：计算区域inum_a的左、右边界的列索引
		    if(ib(inum_a,i) /= 3) then      ! 非interior边界
		        do jpre = gcb(inum_a,3)+1,gcb(inum_a,4)-1       ! jpre：ipre边界上的节点（不包括边界端点）
		            ro_boundary(ipre,jpre) = 2.
		        end do
		    end if
		end do
		
		! 获得每个子计算区域的四个边界交点处的ro_boundary
		do i = 1,2
		    ipre = gcb(inum_a,i)
		    do j = 1,2
		        jpre = gcb(inum_a,j+2)
		        ro_boundary(ipre,jpre) = ro_boundary(ipre,jpre)+1.
		    end do
		end do
	end do
	
	do i=0,nzt
	    do j=0,nrt
	        if(ro_boundary(i,j)>0) ro_boundary(i,j)=4./ro_boundary(i,j)
	    end do
	end do
    
    ! cmprg_gradient：计算边界电场时的判断条件，同时作为计算Poisson方程的判断条件之一，不参与实际计算。
    ! 对于任意子计算区域的非interior类型边界，cmprg_gradient = 2, 4, 5, 6 for left, right, lower, upper boundaries, respectively.
    ! 对于interior类型边界节点和内节点，cmprg_gradient = 1.
    ! 注：任一非interior类型边界角点可以同时属于2个边界，因此cmprg_gradient要叠加。
    !     例如，角点即属于左边界，又是属于下边界，则cmprg_gradient = 2+5 = 7。
	cmprg_gradient = 0
	do i = 1,num_a
		do j = 1,4
			if(ib(i,j) /= 3) then       ! 非interior边界
				if(j < 3) then        ! 子计算区域的左、右边界
					do k = 0,nr(i)
						ipre = gcb(i,j)
						jpre = gcb(i,3)+k
						if (ro_boundary(ipre,jpre) == 2) then        ! 非边界交点的网格节点
							if(j == 1) cmprg_gradient(ipre,jpre) = 2
							if(j == 2) cmprg_gradient(ipre,jpre) = 4
						else if (ro_boundary(ipre,jpre) /= 1) then      ! 边界交点
							if(j == 1) cmprg_gradient(ipre,jpre) = cmprg_gradient(ipre,jpre)+2
							if(j == 2) cmprg_gradient(ipre,jpre) = cmprg_gradient(ipre,jpre)+4
						end if
					end do
				else        ! 子计算区域的上、下边界
					do k = 0,nz(i)
						ipre = gcb(i,1)+k
						jpre = gcb(i,j)
						if(ro_boundary(ipre,jpre) == 2) then        ! 非边界交点的网格节点
							if(j == 3) cmprg_gradient(ipre,jpre) = 5
							if(j == 4) cmprg_gradient(ipre,jpre) = 6
						else if (ro_boundary(ipre,jpre) /= 1) then      ! 边界交点
							if(j == 3) cmprg_gradient(ipre,jpre) = cmprg_gradient(ipre,jpre)+5
							if(j == 4) cmprg_gradient(ipre,jpre) = cmprg_gradient(ipre,jpre)+6
						end if
					end do
				end if
			end if
		end do
		
	end do
	
	! +++++++++++2009-05-08+++++++++++++++++++++
	do i = 1,num_a
	    do j = 0,nz(i)        ! 子计算区域的横向内节点
	        do k = 0,nr(i)        ! 子计算区域的径向内节点
				ipre = gcb(i,1)+j
				jpre = gcb(i,3)+k
				if (cmprg_gradient(ipre,jpre) == 0) cmprg_gradient(ipre,jpre) = 1
			end do
		end do
    end do
    ! +++++++++++END++++++++++++++++++++++++++++
    
    ! +++++++++++2009-05-15+++++++++++++++++++++
    ! cmprg_bdtype：计算Poisson方程的判断条件之一，标识了节点所在的边界类型，不参与实际计算。
    ! 对于任意子计算区域的非interior类型边界，cmprg_bdtype = ib(i,j)。
    ! 对于interior类型边界和内节点，cmprg_bdtype = 3。
    ! 注：任一非interior类型边界的端点可以同时属于2个子边界，因此cmprg_bdtype要叠加。
    ! 叠加例子：若节点所在的一边属于dielectric surface边界，另一边属于periodic边界，则cmprg_bdtype = 24 or 42，即"2"与"4"的组合。
	cmprg_bdtype = 0
	do i = 1,num_a
		do j = 1,4
			if(ib(i,j) /= 3) then       ! 非interior边界
			    if(ib(i,j)==2 .OR. ib(i,j)==4 .OR. ib(i,j)==5) then     ! 目前只考虑了除了interior边界之外的2、4、5类边界条件
			        if(j < 3) then        ! 子计算区域的左、右边界
					    do k = 0,nr(i)
						    ipre = gcb(i,j)
						    jpre = gcb(i,3)+k
						    if (ro_boundary(ipre,jpre) == 2) then        ! 非边界交点的网格节点
						        cmprg_bdtype(ipre,jpre) = ib(i,j)
						    else if (ro_boundary(ipre,jpre) /= 1) then      ! 边界交点
						        if(cmprg_bdtype(ipre,jpre) == 0) then
						            cmprg_bdtype(ipre,jpre) = cmprg_bdtype(ipre,jpre) +ib(i,j)*10
						        else
						            if (ib(i,j)<10) then
						                cmprg_bdtype(ipre,jpre) = cmprg_bdtype(ipre,jpre) +ib(i,j)
						            else
						                cmprg_bdtype(ipre,jpre) = cmprg_bdtype(ipre,jpre)*10 +ib(i,j)
						            end if
						        end if
						    end if
						end do
				    else        ! 子计算区域的上、下边界
					    do k = 0,nz(i)
						    ipre = gcb(i,1)+k
						    jpre = gcb(i,j)
						    if(ro_boundary(ipre,jpre) == 2) then        ! 非边界交点的网格节点
							    cmprg_bdtype(ipre,jpre) = ib(i,j)
						    else if (ro_boundary(ipre,jpre) /= 1) then      ! 边界交点
						        if(cmprg_bdtype(ipre,jpre) == 0) then
						            cmprg_bdtype(ipre,jpre) = cmprg_bdtype(ipre,jpre) +ib(i,j)*10
						        else 
						            if (ib(i,j)<10) then
						                cmprg_bdtype(ipre,jpre) = cmprg_bdtype(ipre,jpre) +ib(i,j)
						            else
						                cmprg_bdtype(ipre,jpre) = cmprg_bdtype(ipre,jpre)*10 +ib(i,j)
						            end if
						        end if
						    end if
						end do
				    end if
				else
				    write (*, "('Unhandled Boundary Condition ', I8)") ib(i,j)
				    stop        ! can comment if not necessary
				end if
			end if
		end do
	end do
	
	do i = 1,num_a
	    do j = 0,nz(i)        ! 子计算区域的横向内节点
	        do k = 0,nr(i)        ! 子计算区域的径向内节点
				ipre = gcb(i,1)+j
				jpre = gcb(i,3)+k
				if (cmprg_bdtype(ipre,jpre) == 0) cmprg_bdtype(ipre,jpre) = 3
			end do
		end do
    end do
    
    ro_boundary_flag = 0	
    do i=0,nzt
        do j=0,nrt
            if(cmprg_gradient(i,j) == 1) ro_boundary(i,j) = 1.
            if ( abs(ro_boundary(i,j)-1.0)<=float_error ) then
                ro_boundary_flag(i,j) = 1	! 内点
            else if ( abs(ro_boundary(i,j)-2.0)<=float_error ) then
                ro_boundary_flag(i,j) = 2	! 边界点
            else if ( abs(ro_boundary(i,j)-4./3.)<=float_error ) then
                ro_boundary_flag(i,j) = 3	! 外角点
            else if ( abs(ro_boundary(i,j)-4.0)<=float_error ) then
                ro_boundary_flag(i,j) = 4	! 内角点
            end if
        end do
    end do
    ! +++++++++++END++++++++++++++++++++++++++++    	

	x_solver_flag = 0
	do j = 0,nrt
		if (ro_boundary_flag(0,j)==0) then
			x_solver_flag(j,1) = 0
		else
			if (cmprg_bdtype(0,j)==4) then
				x_solver_flag(j,1) = 1
			else if (cmprg_bdtype(0,j)==24 .OR. cmprg_bdtype(0,j)==42) then
				if (cmprg_bdtype(1,j)==4) then
					x_solver_flag(j,1) = 0
				else
					x_solver_flag(j,1) = 1
				end if
			else
				x_solver_flag(j,1) = 0
			end if
		end if
	end do
	
	y_solver_flag = 0
	do i = 0,nzt
		if (ro_boundary_flag(i,0)==0) then
			y_solver_flag(i,1) = 0
		else
			if (cmprg_bdtype(i,0)==4) then
				y_solver_flag(i,1) = 1
			else if (cmprg_bdtype(i,0)==24 .OR. cmprg_bdtype(i,0)==42) then
				if (cmprg_bdtype(i,1)==4) then
					y_solver_flag(i,1) = 0
				else
					y_solver_flag(i,1) = 1
				end if
			else
				y_solver_flag(i,1) = 0
			end if
		end if
	end do

!	call viewer(ro_boundary,0,imax1,0,imax2,"ro_boundary")
!	call viewer(real(cmprg_gradient),0,imax1,0,imax2,"cmprg_gradient")
!	call viewer(real(cmprg_bdtype),0,imax1,0,imax2,"cmprg_bdtype")
!	call viewer(real(x_solver_flag),0,imax,1,1,"x_solver_flag")
!	call viewer(real(y_solver_flag),0,imax,1,1,"y_solver_flag")
	
    ! read magnetic field
	allocate(B((nzt+1)*(nrt+1), 2))
	B = 0.
	
	open(unit=10, file='.\input\B.dat')
!	READ(10,*) B
	close(10)

    bzt = 0.
    brt = 0.
    bt = 0.
    
    k = 0
    DO i = 0, nzt
		DO j = 0, nrt
			k = k+1
			if(cmprg_gradient(i,j) > 0) then
!			    bzt(i,j) = B(k,1)
!			    brt(i,j) = B(k,2)
			    bzt(i,j) = 0.
			    brt(i,j) = 0.02
			    bt(i,j) = sqrt(bzt(i,j)*bzt(i,j)+brt(i,j)*brt(i,j))
			end if
		END DO
	END DO
	
	deallocate(B)
    
	bmax = MAXVAL(ABS(bt))
	
    ! if igas = 2, the cross section of krypton should be imported
    if(igas == 2) then
    	open(unit=10, file='.\input\krypton.dat')
		    READ(10,*) krypton
	    close(10)
	    sigmaElastic_kr = krypton(:,2)
	    sigmaExc_kr = krypton(:,3)
	    sigmaIz_kr = krypton(:,4)
    end if
   
    ! initialize electric field
	ez = 0.
	er = 0.
	phi = 0.
	phi0 = 0.

    LMDD = SQRT(EPSILON0*kb*11600*t_int(1)/(e*e*nero)) ! debye length
	OMIGP = SQRT(E*NERO*E/EPSILON0/M(1)) ! plasma frequency
	
	! 确定泊松方程被改变的倍数，德拜长度和等离子震荡频率相应地被改变
	if (dzr > lmdd) then
	    pfactor = (dzr/lmdd)**2
		WRITE(*,"('pfactor = ', F10.4)") pfactor
	    LMDD = LMDD*sqrt(pfactor)
	    OMIGP = OMIGP/sqrt(pfactor)
	end if
	
	VETH = SQRT(kb*T_int(1)*11600/M(1)) ! thermal velocity

	omigc(1:nsp) = e*bmax/(m(1:nsp)) ! cycle frequency
	write(*,*) "omigc = ", omigc(1:nsp)
    dt = 0.35*omigp/omigc(1)
    if(dt > 0.3) dt = 0.3
	dt = 0.1
    write(*,*) "dt = ", dt
	octp = omigc(1)/omigp ! 电子回旋频率和电子震荡频率的比，用于电子运动的计算

	t_parameter = 0.5*m*veth**2/e

	! parameter of ion
	vir_ion0 = SQRT(kb*t_int(1)*11600/m(2)) ! bohm velocity
	fluxqi = e*nero*vir_ion0*dt/OMIGP ! 一个时间步长内离子的电荷个数
	
	! Normalization
	cb = cb/lmdd
	lcb = lcb/lmdd
	z1 = z1/lmdd
	z2 = z2/lmdd
	r1 = r1/lmdd
	r2 = r2/lmdd
	r0 = r0/lmdd
	dz = dz/lmdd
	dr = dr/lmdd
	dzr = dzr/lmdd
	drr = drr/lmdd

	if(bmax > 0 ) then
		brt = brt/bmax
		bzt = bzt/bmax
	else
		brt = 0.
		bzt = 0.
	end if
    
	ez = ez*lmdd/t_int(1)
	er = er*lmdd/t_int(1)
	phi = phi/t_int(1)
	phi0 = phi0/t_int(1)
	
	DO j = 0, nrt
		! 依赖于坐标系
		volume_cell(j) = dzr *drr *lmdd**2
	END DO
	
	! now initialize paticle number, location and velocity
	z = 0.
    r = 0.
	r_ini = 0.
    vz = 0.
    vr = 0.
    vtheta = 0.
	t = 0.
	collectq = 0.

	SEE_flag = 0

	! 依赖于坐标系
	volumet(:) = (cb(:,4)-cb(:,3)) *(cb(:,2)-cb(:,1)) *lmdd**2
!	do insp = 1, nsp
!		num(insp) = 5.e6 ! 假象系统稳定时每种粒子的真实个数
!		! 计算每个巨粒子代表多少个真实粒子
!		do i = 1, num_a
!			fnum(i) = nero*sum(volumet)/num(insp)
!		end do
!	end do

	fnum = 5.e6

	num = 0
	
	continue

END SUBROUTINE INITIAL

subroutine loadv(v, energy, insp)
	use constant
	use global
	real pacc,v,belta,energy
	integer insp

	belta = sqrt(m(insp)/(2*kb*energy*11600))

!	if(insp == 1) then
!		belta = sqrt(m(1)/(2*kb*energy*11600))
!	else
!		belta = sqrt(m(2)/(2*kb*energy*11600))
!	end if

100	call random(ranum)

	v = (6*ranum-3)/belta ! Maximum velocity values injected are six times the thermal velocity
	pacc = exp(-v**2*belta**2)
	call random(ranum1)
	if (ranum1 > pacc) goto 100

end subroutine loadv
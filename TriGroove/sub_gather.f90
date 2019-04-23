SUBROUTINE gather
	USE constant
	use global
	use wall
	use ion
	use diag
	use neutral
	use grid_para
	use mod_dadi

	implicit none

	real ra,rb,rc,area,para,p1,p2,p3,p4
	real  dt1,dt2,dtmin,DT_wall(2)
	integer INUM_A,ipret,jpret,i,j,k,insp,ICOUNT,IBOUND,inspa,inspi

	do insp = 1, nsp
		if(insp==1 .OR. insp==2) then
			! 统计宏观参数前需要清零
			rho(insp,:,:) = 0.

			do j = 1, num(insp)
				INUM_A = 0
				DO i = 1, NUM_A
					IF(Z(insp,j) > CB(i,1) .AND. Z(insp,j) < CB(i,2) .AND.  &
						R(insp,j) > CB(i,3) .AND. R(insp,j) < CB(i,4)) INUM_A = i
				END DO
				
				IF(INUM_A == 0)  then
					z(insp,j) = -100 ! 如果粒子超出计算区域（这是不合理的），则将该粒子做上标记，循环完后统一注销
					goto 100
				end if

				ipret = INT((z(insp,j) -z1) /dzr)
				jpret = INT((r(insp,j) -r1) /drr)
				if(ipret>=0 .AND. ipret<=nzt .AND. jpret>=0 .AND. jpret<=nrt) then
					IF(ipret == nzt) then
						ipret = ipret-1
						z(insp,j) = cb(inum_a, 2) - dzr*0.001
					end if
					IF(jpret == nrt) then
						jpret = jpret-1
						r(insp,j) = cb(inum_a, 4) - drr*0.001
					end if
					call co_ptogt(insp,j,inum_a,ipret,jpret,p1,p2,p3,p4) ! 计算各种密度的系数p1,p2,p3,p4
					para = 1.
					call ptogt(rho,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
					call ptogt(rho_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
					! rho 与 rho_g1 的区别：
					! rho 是当前时刻的密度
					! rho_g1 是累计了 case_output 步长之后的密度，因此，rho_g1 约等于 rho*case_output。
	
					! 计算网格节点上的通量密度rhovz_g1(轴向)、rhovr_g1(径向)、rhovtheta_g1(周向)以及能量密度rhoegy_g1
					para = vz(insp,j)
					call ptogt(rhovz_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
					para = vr(insp,j)
					call ptogt(rhovr_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
					if (insp == 1) then
				        para = vtheta(j)
				        call ptogt(rhovtheta_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)

						if (SEE_flag(j) == 1) then
							call ptogt(rho_g1_SEE,nsp,imax1,imax2,j,insp,1.,ipret,jpret,p1,p2,p3,p4)
							call ptogt(rhovtheta_g1_SEE,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
						end if

				    end if
	
					para = t(insp,j)
					call ptogt(rhoegy_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
	
					if (insp==1) then
						para = vz(insp,j)**2 +vr(insp,j)**2 +vtheta(j)**2
					else
						para = vz(insp,j)**2 +vr(insp,j)**2
					end if
					call ptogt(rhovt2_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)

					if (insp==1) then
						para = t_parameter(insp) * (vz(insp,j)**2 +vtheta(j)**2)
					else
						para = t_parameter(insp) * vz(insp,j)**2
					end if
					call ptogt(rhoegy_g1_normal,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)
					
					para = sqrt(para)
					call ptogt(rhovt_g1,nsp,imax1,imax2,j,insp,para,ipret,jpret,p1,p2,p3,p4)

				else
					z(insp,j) = -100 ! 如果粒子超出计算区域（这是不合理的），则将该粒子做上标记，循环完后统一注销
	
				end if

100				continue

			end do

			j = 1
			do while (j <= num(insp))
				if(z(insp,j) < 0) call remove(insp,j)
				j = j+1
			end do

		end if
	end do

	! 计算绝缘壁面的表面电荷密度COLLECTQ
	do inum_a = 1, num_a
		do i = 1, 4
			if(ib(inum_a,i) == 2) then
				if(i==1 .OR. i==2) then
					area = dr(inum_a)*lmdd
					do j=0, nr(inum_a)
						if(nsp == 1) then
							COLLECTQ(inum_a,i,j) = COLLECTQ(inum_a,i,j) &
													+ COLLECT(1,inum_a,i,j)*Q(1)/area*lmdd/(T_int(1)*EPSILON0*pfactor) &
													+ FLUXQI(inum_a,j)*lmdd/(T_int(1)*EPSILON0*pfactor)
						else
							inspi = 2
							COLLECTQ(inum_a,i,j) = COLLECTQ(inum_a,i,j) &
													+ (COLLECT(1,inum_a,i,j)*Q(1)+COLLECT(inspi,inum_a,i,j)*Q(inspi))/area &
													*lmdd/(T_int(1)*EPSILON0*pfactor)
						end if
					end do
				end if
				
				if(i==3 .OR. i==4) then
					area = dz(inum_a)*lmdd
					do j = 0,nz(inum_a)
						if(nsp==1) then
							COLLECTQ(inum_a,i,j) = COLLECTQ(inum_a,i,j) &
													+ COLLECT(1,inum_a,i,j)*Q(1)/area*lmdd/(T_int(1)*EPSILON0*pfactor) &
													+ FLUXQI(inum_a,j)*lmdd/(T_int(1)*EPSILON0*pfactor)
						else
							inspi = 2
							COLLECTQ(inum_a,i,j) = COLLECTQ(inum_a,i,j) &
													+ (COLLECT(1,inum_a,i,j)*Q(1)+COLLECT(inspi,inum_a,i,j)*Q(inspi))/area &
													*lmdd/(T_int(1)*EPSILON0*pfactor)
						end if
					end do
				end if
			end if
		end do
	end do

END SUBROUTINE gather

SUBROUTINE REMOVE(insp,j)
	! 粒子消失的几种情况：从出口喷出，电子和离子在壁面的积累，电子和离子打到阳极
	
	use neutral
	use global
	use diag
	implicit none

	integer insp, j

	Z(insp,j) = Z(insp,num(insp))
	R(insp,j) = R(insp,num(insp))
	
	VZ(insp,j) = VZ(insp,num(insp))
	VR(insp,j) = VR(insp,num(insp))
	if(insp == 1) then
		Vtheta(j) = Vtheta(num(insp))
		SEE_flag(j) = SEE_flag(num(insp))
		SEE_flag(num(insp)) = 0
	end if
	
	if(insp == 3) fnuma(j) = fnuma(num(insp))
	
	if(insp == 1) then
		T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2+Vtheta(j)**2)
	else
		T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2)
	end if
	
	num(insp) = num(insp) - 1
	j = j - 1

END SUBROUTINE REMOVE

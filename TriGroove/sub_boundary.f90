subroutine boundary(insp,ibound,j,inum_a)
	use constant
	use global
	use wall
	use neutral
	use diag
	implicit none
	
	integer ibound,i,j,inum_a,inum_a1,k,insp,jpre,inspa,new_particle,ipre
	real tez,ter,vmp,energy,v1,v2,v3

	IF (IBOUND == 1) THEN
		IF (z(insp,j) <= CB(INUM_A,1)) K = 1 ! lower boundary for z direction
		IF (z(insp,j) >= CB(INUM_A,2)) K = 2 ! upper boundary
		
		AT(insp,j) = (z(insp,j)-CB(INUM_A,K))/vz(insp,j) ! 穿越边界之后剩余的运动时间
		
		IF (IB(INUM_A,K)==1 .OR. IB(inum_a,k)==11) THEN ! 真空边界1与导体边界11，粒子消失
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			z(insp,j) = -100
			if(ib(inum_a,k) == 11) then
				if(insp==1 .OR. insp==2) then
					conductor_charge = conductor_charge + fnum(1)*q(insp)
				end if
			end if

		ELSE IF (IB(INUM_A,K) == 2 .OR. ib(inum_a,k)==8) THEN ! 绝缘壁面边界2与绝缘并固定电势壁面边界8
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			z(insp,j) = CB(INUM_A, K) +0.0001*DZ(INUM_A)*(3.-2.*K)
			r(insp,j) = r(insp,j) - vr(insp,j)*AT(insp,j)
			CALL dielectric(INUM_A,j,IBOUND,K,insp)

		ELSE IF (IB(INUM_A,K) == 3)	THEN ! 内部边界3
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			continue
		
		ELSE IF(IB(INUM_A,K) == 4) THEN ! 周期性边界4
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			ICOLLIDE_WALL(insp,j) = 1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if(inum_a==1 .AND. k==1) then
				z(insp,j) = CB(3,2) - 0.0001*DZ(3)
			else if(inum_a==3 .AND. k==2) then
				z(insp,j) = CB(1,1) + 0.0001*DZ(1)
			end if

			! plane
!			if(inum_a==1 .AND. k==1) then
!				z(insp,j) = CB(1,2) - 0.0001*DZ(1) 
!			else if(inum_a==1 .AND. k==2) then
!				z(insp,j) = CB(1,1) + 0.0001*DZ(1)
!			end if
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			r(insp,j) = r(insp,j) - vr(insp,j)*AT(insp,j)
			IF(r(insp,j) < CB(INUM_A,3)) r(insp,j) = CB(INUM_A,3)+0.0001*DR(INUM_A)
			IF(r(insp,j) > CB(INUM_A,4)) r(insp,j) = CB(INUM_A,4)-0.0001*DR(INUM_A)

		ELSE if(ib(inum_a,k) == 5) then ! 入口边界5
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			z(insp,j) = -100

		ELSE IF(ib(inum_a,k)==6 .OR. ib(inum_a,k)==9) then ! 镜面反射边界6以及轴对称边界9
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			z(insp,j) = 2*CB(INUM_A, K)-z(insp,j)
			vz(insp,j) = -vz(insp,j)

		else if(IB(INUM_A,K) == 7)	THEN ! 阳极边界7
			nout(insp,inum_a,k) = nout(insp,inum_a,k) + 1
			z(insp,j) = -100

		else if(ib(inum_a,k)==10) then ! 无穷远边界10
			nout(insp,inum_a,k) = nout(insp,inum_a,k) +1
			
			if(insp==1) then
				! 电子被反射回到计算区域
				ICOLLIDE_WALL(insp,j) = 1
				r(insp,j) = r(insp,j) - vr(insp,j)*AT(insp,j)
				z(insp,j) = CB(INUM_A, K) + 0.0001*DZ(INUM_A)*(3.-2.*K)
				vmp = sqrt(2.*kb*11600*t(insp,j)/m(insp))/veth
				call RVELC(v1, v2, v3, VMP, k, 1)
				vr(insp,j) = v1
				vtheta(j) = v2
				vz(insp,j) = v3
			else
				z(insp,j) = -100
			end if
		end if
	END IF

	IF (IBOUND == 2) THEN
		IF (r(insp,j) <= CB(INUM_A, 3)) K = 1
		IF (r(insp,j) >= CB(INUM_A, 4)) K = 2
		
		AT(insp,j) = (r(insp,j)-CB(INUM_A, K+2))/vr(insp,j)
		
		IF(IB(INUM_A,K+2) == 1 .OR. IB(inum_a,K+2) == 11) THEN  ! 真空边界1与导体边界11，粒子消失
			nout(insp,inum_a,k+2) = nout(insp,inum_a,k+2) + 1
			z(insp,j)=-100
			if(ib(inum_a,k+2) == 11) then
				if(insp==1 .OR. insp==2) then
					conductor_charge = conductor_charge + fnum(1)*q(insp)
				end if
			end if

		ELSE IF(IB(INUM_A,K+2) == 2 .OR. ib(inum_a,k+2) == 8) THEN ! 绝缘壁面边界2与绝缘并固定电势壁面边界8
			nout(insp,inum_a,k+2) = nout(insp,inum_a,k+2) + 1
			r(insp,j) = CB(INUM_A, K+2) + 0.0001*DR(INUM_A)*(3.-2.*K)
			z(insp,j) = z(insp,j) - vz(insp,j)*AT(insp,j)
			CALL dielectric(INUM_A,j,IBOUND,K,insp)

		ELSE IF(IB(INUM_A,K+2) == 3) THEN ! 内部边界3
			nout(insp,inum_a,k+2) = nout(insp,inum_a,k+2) + 1
			continue

		ELSE if(ib(inum_a,k+2) ==4) then ! 周期性边界4
			continue
			
		ELSE IF(IB(INUM_A,K+2) == 5) THEN ! 入口边界5
			nout(insp,inum_a,k+2) = nout(insp,inum_a,k+2) + 1
			z(insp,j) = -100

		ELSE IF(IB(INUM_A,K+2) == 6 .OR. ib(inum_a,k+2)==9) THEN ! 镜面反射边界6以及轴对称边界9
			nout(insp,inum_a,k+2) = nout(insp,inum_a,k+2) + 1
			r(insp,j) = 2*CB(INUM_A, K+2) - r(insp,j)
			vr(insp,j) = -vr(insp,j)

		else if(IB(INUM_A,K+2) == 7) THEN ! 阳极边界7
			nout(insp,inum_a,k+2) = nout(insp,inum_a,k+2) +1
			z(insp,j) = -100

		else if(ib(inum_a,k+2) == 10) then ! 无穷远边界10
			nout(insp,inum_a,K+2) = nout(insp,inum_a,K+2) +1
			
			if(insp==1) then
				! 电子被反射回到计算区域
				ICOLLIDE_WALL(insp,j) = 1
				z(insp,j) = z(insp,j)-vz(insp,j)*AT(insp,j)
				r(insp,j) = CB(INUM_A, K+2)+0.0001*DR(INUM_A)*(3.-2.*K)
				vmp = sqrt(2.*kb*11600*t(insp,j)/m(insp))/veth
				call RVELC(v1,v2,v3,VMP,k,1)
				vz(insp,j) = v1
				vtheta(j) = v2
				vr(insp,j) = v3
			else
				z(insp,j) = -100
			end if

		end if
	END IF

end subroutine boundary

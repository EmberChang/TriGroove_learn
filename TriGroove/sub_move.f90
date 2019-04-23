SUBROUTINE MOVE
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

	! 每次粒子循环前宏观参数清零
	COLLECT = 0 ! particle on wall every time step
	ICOLLIDE_WALL = 1 ! 1 means collide with boundary
	AT = DT

	do insp = 1, nsp
		if(insp==1 .OR. insp==2) then
			j = 1
			DO while (j <= num(insp))
				particle_add_num(insp) = 0
				
				! Determine which area is this particle in
				do while(ICOLLIDE_WALL(insp,j) == 1) ! 碰撞后ICOLLID_WALL的值改变
					ICOLLIDE_WALL(insp,j) = 0
					INUM_A = 0
					DO i = 1, NUM_A
						IF(Z(insp,j) > CB(i,1) .AND. Z(insp,j) < CB(i,2) .AND.  &
							R(insp,j) > CB(i,3) .AND. R(insp,j) < CB(i,4)) INUM_A = i
					END DO
					
					IF(INUM_A == 0)  then
						z(insp,j) = -100 ! 如果粒子超出计算区域（这是不合理的），则将该粒子做上标记，循环完后统一注销
						goto 100 ! remove particle and goto the end
					end if
					
					! 计算电磁力对带电粒子的运动
					call update_v_p(insp,j,inum_a)
					
					! the following deal with boundary condition
					IBOUND = 0 ! denote boundary,1 for z direction ,2 for r direction
					ICOUNT = 0 ! denote collide with boundary
					IF (Z(insp,j) <= CB(INUM_A, 1) .OR. Z(insp,j) >= CB(INUM_A, 2)) THEN
						IBOUND = 1
						ICOUNT = ICOUNT+1
					END IF
					IF (R(insp,j) <= CB(INUM_A, 3) .OR. R(insp,j) >= CB(INUM_A, 4)) THEN
						IBOUND = 2
						ICOUNT = ICOUNT+1
					END IF
					
					IF (ICOUNT > 1) THEN
						! collide with two boundaries
						DO K=1,2
							DT_wall(K)=5.*at(insp,j)
						END DO
						
						IF (ABS(vz(insp,j)) > 1.E-6) THEN
							DT1=(CB(INUM_A, 1)-z(insp,j))/vz(insp,j)
							DT2=(CB(INUM_A, 2)-z(insp,j))/vz(insp,j)
							DT_wall(1)=DT1
							IF (DT2 > DT1) DT_wall(1)=DT2
						END IF
						
						IF (ABS(vr(insp,j)) > 1.E-6) THEN
							DT1=(CB(INUM_A, 3)-r(insp,j))/vr(insp,j)
							DT2=(CB(INUM_A, 4)-r(insp,j))/vr(insp,j)
							DT_wall(2)=DT1
							IF (DT2 > DT1) DT_wall(2)=DT2
						END IF

						! DTMIN is negative, 
						! its absolute value represents the rest time 
						! that the particle should move for after colliding with the wall
						DTMIN=5.*DT
						DO K=1,2
							IF (DT_wall(K) < DTMIN) THEN
								DTMIN=DT_wall(K)
								IBOUND=K ! IBOUND now indicates the first side that is crossed
							END IF
						END DO
						
						ICOUNT=1 ! now only deal with colliding with one boundary
					END IF
					
					! only collide with one boundary
					IF(ICOUNT == 1)   call boundary(insp,ibound,j,inum_a)
				END DO

100				continue
				
				num(insp) = num(insp) + particle_add_num(insp) ! 将新增加的粒子计入总数
				
				j = j + 1
			END DO
		end if
	end do
	
END SUBROUTINE MOVE

subroutine update_v_p(insp,j,inum_a)
	use global
	use magnetic
	use diag
	implicit none

	integer insp,i,j,inum_a,ipret,jpret
	character aa,bb
	real F,octp_dt,Q_DT,ezpre,erpre,bzpre,brpre
	real A(3),Tm(3),S(3),uprime(3),a1,a2,a3,a4
	real eznew,ernew,EthetaPRE
	real Rbro,thetabro,rpha,rbro1,cosrpha,sinrpha,vrbro,vthetabro

	EthetaPRE = 2.e4*lmdd/t_int(1)
!	...... or ......
!	EthetaPRE = 0.

	if(insp==1 .OR. insp==2) then ! 计算电磁力对带电离子的作用
		F = half_q_m(insp) * AT(insp,j)
		Q_DT = Q1(insp) / AT(insp,j)
		ipret = int((z(insp,j) - z1) / dzr)
		jpret = int((r(insp,j) - r1) / drr)
		if(ipret == nzt) then
			ipret = ipret - 1
			z(insp,j) = z(insp,j) - dzr*0.0001
		end if
		if(jpret == nrt) then
			jpret = jpret - 1
			r(insp,j) = r(insp,j) - drr*0.0001
		end if
		call co_gtop(insp, j, ipret, jpret, a1, a2, a3, a4)
		call gtop(insp, j, ezpre, ez, imax1, imax2, ipret, jpret, a1, a2, a3, a4)
		call gtop(insp, j, erpre, er, imax1, imax2, ipret, jpret, a1, a2, a3, a4)

		if(insp == 1) then
			octp_dt = octp*at(insp,j)
			call gtop(insp, j, bzpre, bzt, imax1, imax2, ipret, jpret, a1, a2, a3, a4)
			call gtop(insp, j, brpre, brt, imax1, imax2, ipret, jpret, a1, a2, a3, a4)
			Tm(1) = octp_dt * bzpre/2
			Tm(2) = octp_dt * brpre/2
			Tm(3) = 0
			S(1) = 2*Tm(1)/(1+Tm(1)*Tm(1))
			S(2) = 2*Tm(2)/(1+Tm(2)*Tm(2))
			S(3) = 2*Tm(3)/(1+Tm(3)*Tm(3))
		end if
		
		A(1) = F*EZPRE
		A(2) = F*ERPRE
		A(3) = F*EthetaPRE

		! half acceleration of axial and radial velocity
		VZ(insp,j) = VZ(insp,j) +A(1)
		VR(insp,j) = VR(insp,j) +A(2)
		
		if(insp==1) then ! 电子在磁场作用下的回旋运动
			vtheta(j) = vtheta(j) +A(3) ! half acceleration of azimuthal velocity

			! rotation
			UPRIME(1) = VZ(insp,j) +VR(insp,j)*Tm(3) -Vtheta(j)*Tm(2) ! x2*v.x3 - x3*v.x2
			UPRIME(2) = VR(insp,j) +Vtheta(j)*Tm(1) -VZ(insp,j)*Tm(3) ! x3*v.x1 - x1*v.x3
			UPRIME(3) = Vtheta(j) +VZ(insp,j)*Tm(2) -VR(insp,j)*Tm(1) ! x1*v.x2 - x2*v.x1
			vz(insp,j) = vz(insp,j) +UPRIME(2)*S(3) -UPRIME(3)*S(2)
			vr(insp,j) = vr(insp,j) +UPRIME(3)*S(1) -UPRIME(1)*S(3)
			vtheta(j) = vtheta(j) +UPRIME(1)*S(2) -UPRIME(2)*S(1)
			
			vtheta(j) = vtheta(j) +A(3) ! half acceleration of azimuthal velocity
		end if

		! half acceleration of axial and radial velocity
		vz(insp,j) = vz(insp,j) +A(1)
		vr(insp,j) = vr(insp,j) +A(2)
	end if

	! update location
	Z(insp,j) = Z(insp,j) +VZ(insp,j)*at(insp,j)
	r(insp,j) = r(insp,j) +vr(insp,j)*at(insp,j)
	
	if(insp == 1) then ! 只考虑电子的径向效应
		T(insp,j) = t_parameter(insp)*(VZ(insp,j)*VZ(insp,j)+VR(insp,j)*VR(insp,j)+Vtheta(j)*Vtheta(j))
	else
		T(insp,j) = t_parameter(insp)*(VZ(insp,j)*VZ(insp,j)+VR(insp,j)*VR(insp,j))
	end if

end subroutine update_v_p

SUBROUTINE RVELC(U,V,W,VMP,K,flag)
! generates two random velocity components U an V in an equilibrium gas with most probable speed VMP
! 半麦克斯韦分布

	USE constant
!	use global
	implicit none

	real A,B,COSPHI,SINPHI,U,V,W,VMP,ranum
	INTEGER K,flag

	if(flag==0) then
		CALL RANDOM(RANUM)
		A = SQRT(-LOG(RANUM))*VMP
		CALL RANDOM(RANUM)
		B = TWOPI*RANUM

		CALL RANDOM(RANUM)
		COSPHI = 1.-RANUM
		SINPHI = SQRT(1-COSPHI**2)
		U = A*SIN(B)*SINPHI
		V = A*COS(B)*SINPHI
		W = (3.-2.*K)*ABS(A*COSPHI)
	else
		A = VMP
		CALL RANDOM(RANUM)
		B = TWOPI*RANUM
		CALL RANDOM(RANUM)
		COSPHI = 1.-RANUM
		SINPHI = SQRT(1-COSPHI**2)
		U = A*SIN(B)*SINPHI
		V = A*COS(B)*SINPHI
		W = (3.-2.*K)*ABS(A*COSPHI)
    END IF

END SUBROUTINE RVELC
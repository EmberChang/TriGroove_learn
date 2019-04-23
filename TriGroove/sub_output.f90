subroutine output

	use constant
	use global
	use magnetic
	use wall
	use diag
	use neutral
	use mcc_para
	use mod_dadi
	use grid_para
	use ion
	use grid
	implicit none

	real Itheta, jetheta_g(0:imax1,0:imax2)
	real Itheta_SEE, jetheta_g_SEE(0:imax1,0:imax2)
	real neo(0:nzt,0:nrt), nio(0:nzt,0:nrt), rho_atom(0:imax1,0:imax2)
	real neoo(0:nzt,0:nrt), nioo(0:nzt,0:nrt)
	real vthetao(0:nzt,0:nrt)
	real vthetao_SEE(0:nzt,0:nrt)
	real vi_abs, angle_cos
	real egy_total(nsp), ave_egy(nsp)
	real SEE_Inc(0:imax1,0:imax2), SEE_Emit(0:imax1,0:imax2), SEE_Coef(0:imax1,0:imax2)
	character(len=80) :: filename,form
	integer i,j,k,insp,inum_a,inspa,ipre,jpre,ipret,jpret,mc

	if(iloop==1 .AND. state==1) then
        open(unit=40,file='.\output\FIELD1.dat')
            WRITE(40,1) ILOOP
1	        FORMAT('TITLE = "STEP =',I6,'"')
            WRITE(40,"(A150)")'VARIABLES = "Z" "R" "Na" '
        close(40)
    else if(state==2 .AND. ifield2==0) then
        ifield2=ifield2+1
        open(unit=40,file='.\output\FIELD201.dat')
            WRITE(40,201) ILOOP
201	        FORMAT('TITLE = "STEP =',I6,'"')
            WRITE(40,"(A150)")' VARIABLES = "Z" "R" "Ne" "Ne_SEE" "Ni" '
        close(40)

		open(unit=40,file='.\output\FIELD202.dat')
            WRITE(40,202) ILOOP
202	        FORMAT('TITLE = "STEP =',I6,'"')
            WRITE(40,"(A150)")' VARIABLES = "Z" "R" "vez" "ver" "vetheta" "jtheta" "jtheta_SEE" '
        close(40)

		open(unit=40,file='.\output\FIELD203.dat')
            WRITE(40,203) ILOOP
203	        FORMAT('TITLE = "STEP =',I6,'"')
            WRITE(40,"(A150)")' VARIABLES = "Z" "R" "phi" "Te" "egy_e" "egy_i" '
        close(40)

	end if

	ez_g1 = ez_g1 + ez
	er_g1 = er_g1 + er
	phi_g1 = phi_g1 + phi

	egy_total = 0.
	ave_egy = 0.
	do i = 1, nsp
		do j = 1, num(i)
			egy_total(i) = egy_total(i) + t(i,j)
		end do
		if(num(i)>0) ave_egy(i) = egy_total(i) / num(i)
	end do

	if(mod(iloop,case_output)==0) then ! 每隔case_output步保存计算结果
		WRITE (*,*) 'writing restart and output files'
		if(state == 1) then
			inspa = 3
			ifield1 = ifield1+1
			OPEN(4, FILE='.\output\case.res', FORM='UNFORMATTED')
				WRITE(4) iloop, num(inspa), z(inspa,1:num(inspa)), r(inspa,1:num(inspa)), & 
							r_ini(inspa,1:num(inspa)), vz(inspa,1:num(inspa)), vr(inspa,1:num(inspa)), ifield1
			CLOSE(4)
			
		else if(state==2) then
			ifield2 = 1
			write(form,'(i1)') ifield2
			write(filename, *) ".\output\case_1", trim(form), ".res"
			open(4, file=filename, FORM='UNFORMATTED')
            WRITE(4) iloop, num, ((z(i,j),j=1,num(i)),i=1,nsp), ((r(i,j),j=1,num(i)),i=1,nsp), &
						((r_ini(i,j),j=1,num(i)),i=1,nsp), ((vz(i,j),j=1,num(i)),i=1,nsp), &
						((vr(i,j),j=1,num(i)),i=1,nsp), vtheta(1:num(1)), ((t(i,j),j=1,num(i)),i=1,nsp), &
                        fnuma(1:num(nsp)), fnum, collectq, rho, conductor_charge, ifield2, pfactor, &
                        dt, lmdd, omigp, nout, ez, er
			CLOSE (4)
		end if

		! 速度、能量、温度等参数的统计平均值
		rho_g = 0.
		vz_g = 0.
		vr_g = 0.
		vtheta_g = 0.
		egy_g = 0.
		t_g = 0.
		vt_g = 0.
		vt2_g = 0.
		rho_g_SEE = 0.
		vtheta_g_SEE = 0.
		egy_g_normal = 0.
		do insp = 1, nsp
			do i = 0, nzt
				do j = 0, nrt
					if(rho_g1(insp,i,j) > 0.) then
						rho_g(insp,i,j) = rho_g1(insp,i,j)/case_output
						vz_g(insp,i,j) = rhovz_g1(insp,i,j)*veth*sqrt(mfactor(insp))/rho_g1(insp,i,j)
						vr_g(insp,i,j) = rhovr_g1(insp,i,j)*veth*sqrt(mfactor(insp))/rho_g1(insp,i,j)
						if (insp==1) then
							vtheta_g(insp,i,j) = rhovtheta_g1(insp,i,j)*veth*sqrt(mfactor(insp))/rho_g1(insp,i,j)

							if(rho_g1_SEE(insp,i,j) > 0.) then
								rho_g_SEE(insp,i,j) = rho_g1_SEE(insp,i,j)/case_output
								vtheta_g_SEE(insp,i,j) = rhovtheta_g1_SEE(insp,i,j)*veth*sqrt(mfactor(insp))/rho_g1_SEE(insp,i,j)
							end if
						end if
						
						egy_g(insp,i,j) = rhoegy_g1(insp,i,j)/rho_g1(insp,i,j)
						egy_g_normal(insp,i,j) = rhoegy_g1_normal(insp,i,j)/rho_g1(insp,i,j)

						vt_g(insp,i,j) = rhovt_g1(insp,i,j)*veth*sqrt(mfactor(insp))/rho_g1(insp,i,j)
						vt2_g(insp,i,j) = rhovt2_g1(insp,i,j)*veth**2*mfactor(insp)/rho_g1(insp,i,j)
						t_g(insp,i,j) = m(insp)/(3*e)*(vt2_g(insp,i,j)-vt_g(insp,i,j)**2)
						if (t_g(insp,i,j) < 0) t_g(insp,i,j) = 0
					end if
				end do
			end do
		end do
		ez_g = ez_g1/case_output
		er_g = er_g1/case_output
		phi_g = phi_g1/case_output

		neoo = rho_g(1,:,:)
		nioo = rho_g(2,:,:)
!		call viewer(neoo,0,nzt,0,nrt,"ne_ave")
!		call viewer(nioo,0,nzt,0,nrt,"ni_ave")
!		call viewer(phi_g,0,nzt,0,nrt,"phi_ave")

		! 传导电流及密度的计算
		Itheta = 0.
		jetheta_g = 0.
		Itheta_SEE = 0.
		jetheta_g_SEE = 0.
		do i = 0, nzt
			do j = 0, nrt
				jetheta_g(i,j) = vtheta_g(1,i,j) *rho_g(1,i,j) *e
				Itheta = Itheta +jetheta_g(i,j) *dzr *drr *lmdd**2

				jetheta_g_SEE(i,j) = vtheta_g_SEE(1,i,j) *rho_g_SEE(1,i,j) *e
				Itheta_SEE = Itheta_SEE +jetheta_g_SEE(i,j) *dzr *drr *lmdd**2
			end do
		end do

		vthetao = vtheta_g(1,:,:)
		vthetao_SEE = vtheta_g_SEE(1,:,:)
!		call viewer(vthetao,0,nzt,0,nrt,"vtheta_ave")
!		call viewer(jetheta_g,0,nzt,0,nrt,"je_ave")
!		call viewer(vthetao_SEE,0,nzt,0,nrt,"vtheta_SEE_ave")
!		call viewer(jetheta_g_SEE,0,nzt,0,nrt,"je_SEE_ave")

		! 二次电子发射系数的统计值
		SEE_Inc = 0.
		SEE_Emit = 0.
		SEE_Coef = 0.
		do i = 1, num_a
			do j = 1, 4
				if(ib(i,j) == 2) then
					if(j < 3) then
						do k = 0, nr(i)
							ipre = gcb(i,j)
							jpre = gcb(i,3) +k
							SEE_Inc(ipre,jpre) = SEE_Inc(ipre,jpre) +IncNum(1,i,j,k)
							SEE_Emit(ipre,jpre) = SEE_Emit(ipre,jpre) +EmitNum(1,i,j,k)
						end do
					else
						do k = 0,nz(i)
							ipre = gcb(i,1) +k
							jpre = gcb(i,j)
							SEE_Inc(ipre,jpre) = SEE_Inc(ipre,jpre) +IncNum(1,i,j,k)
							SEE_Emit(ipre,jpre) = SEE_Emit(ipre,jpre) +EmitNum(1,i,j,k)
						end do
					end if
				end if
			end do
		end do

		do i = 0, nzt
			do j = 0, nrt
				if(SEE_Inc(i,j) > 0.) then
					SEE_Coef(i,j) = SEE_Emit(i,j)/SEE_Inc(i,j)
				end if
			end do
		end do

!		call viewer(SEE_Coef,0,nzt,0,nrt,"SEE_Coef")

		! 原子密度
		if(nsp >= 3) then
			inspa=3
			rho_atom = rho_g(inspa,:,:)
		end if


		! 将数据存储到文件		
		if(state==1 .AND. animation==1) then
			open(unit=40,file='.\output\FIELD1.dat',position='append')
				WRITE(40,5) iloop*dt*atomt*neutraltime/omigp,imax1+1,imax2+1
5				FORMAT('ZONE T="Temp. distribution,Time=',E10.5,'s" ','I=',I3,',','J=',I3,',','F=POINT')
				DO j = 0, imax2
					DO k = 0, imax1
						WRITE(40,"(2(F15.10,1X),1(E15.5,1x))") (K*DZr+z1)*lmdd,(j*DRr+r1+r0)*lmdd,rho_atom(k,j)
					END DO
				END DO
				
				write(40,6) (iloop-iloopa)*dt*atomt*neutraltime/omigp,ifield1
6				format('TEXT X=70, Y=90,T="Time=',E10.5,'s",F=COURIER, CS=FRAME, H=2, ZN=',I4)
			CLOSE(40)
			
		else if(state==2 .AND. animation==1) then
			if(iloop>2e5) then
			open(unit=40,file='.\output\FIELD201.dat',position='append')
				WRITE(40,701) iloop*dt/omigp,imax1+1,imax2+1
701				FORMAT('ZONE T="Temp. distribution,Time=',E10.5,'s" ','I=',I3,',','J=',I3,',','F=POINT')
				DO j=0,imax2
					DO k=0,imax1
						WRITE(40,"(2(F15.10,1X), 3(E15.5,1x))") &
							(K*DZr+z1)*lmdd, (j*DRr+r1+r0)*lmdd, rho_g(1,k,j), rho_g_SEE(1,k,j), rho_g(2,k,j)
					END DO
				END DO
				
				write(40,801) iloop*dt/omigp,ifield2
801				format('TEXT X=70, Y=90,T="Time=',E10.5,'s",F=COURIER, CS=FRAME, H=2, ZN=',I4)
			CLOSE(40)

			open(unit=40,file='.\output\FIELD202.dat',position='append')
				WRITE(40,702) iloop*dt/omigp,imax1+1,imax2+1
702				FORMAT('ZONE T="Temp. distribution,Time=',E10.5,'s" ','I=',I3,',','J=',I3,',','F=POINT')
				DO j=0,imax2
					DO k=0,imax1
						WRITE(40,"(2(F15.10,1X), 5(E15.5,1X))") &
							(K*DZr+z1)*lmdd, (j*DRr+r1+r0)*lmdd, vz_g(1,k,j), vr_g(1,k,j), vtheta_g(1,k,j), &
							jetheta_g(k,j), jetheta_g_SEE(k,j)
					END DO
				END DO
				
				write(40,802) iloop*dt/omigp,ifield2
802				format('TEXT X=70, Y=90,T="Time=',E10.5,'s",F=COURIER, CS=FRAME, H=2, ZN=',I4)
			CLOSE(40)

			open(unit=40,file='.\output\FIELD203.dat',position='append')
				WRITE(40,703) iloop*dt/omigp,imax1+1,imax2+1
703				FORMAT('ZONE T="Temp. distribution,Time=',E10.5,'s" ','I=',I3,',','J=',I3,',','F=POINT')
				DO j=0,imax2
					DO k=0,imax1
						WRITE(40,"(2(F15.10,1X), 4(F20.10,1X))") &
							(K*DZr+z1)*lmdd, (j*DRr+r1+r0)*lmdd, &
							phi_g(k,j)*t_int(1), t_g(1,k,j), egy_g(1,k,j), egy_g(2,k,j)
					END DO
				END DO
				
				write(40,803) iloop*dt/omigp,ifield2
803				format('TEXT X=70, Y=90,T="Time=',E10.5,'s",F=COURIER, CS=FRAME, H=2, ZN=',I4)
			CLOSE(40)
			end if	
		end if
				
		open(unit=40,file='.\output\FIELD.dat')
			WRITE(40,"(A400)")'VARIABLES = "Z" "R" "I" "J" "Bz" "Br" & 
											"ne_ave" "ni_ave" "phi_ave" &
											"Te_ave" "egye_ave" "egyi_ave" &
											"egye_n_ave" "egyi_n_ave" "SEE_Coef" & 
											"ez_ave" "er_ave" "viz_ave" "vir_ave" &
											"vez_ave" "ver_ave" "vtheta_ave" "jetheta_g" &
											"ne_SEE_ave" "jetheta_g_SEE" '
			WRITE(40,9) iloop*dt/omigp,imax1+1,imax2+1
9			FORMAT('ZONE T="Temp. distribution,Time=',E10.5,'s" ','I=',I3,',','J=',I3,',','F=POINT')
			DO j=0,imax2
				DO k=0,imax1
					WRITE(40,"(2(F15.10,1X),2(I4,1X),2(F15.10,1X),2(E15.5,1x),7(F15.10,1X),10(E15.5,1x))") &
							(K*DZr+z1)*lmdd, (j*DRr+r1+r0)*lmdd, K, J, bzt(k,j)*bmax, brt(k,j)*bmax, & 
							rho_g(1,k,j), rho_g(2,k,j), phi_g(k,j)*t_int(1), &
							t_g(1,k,j), egy_g(1,k,j), egy_g(2,k,j), &
							egy_g_normal(1,k,j), egy_g_normal(2,k,j), SEE_Coef(k,j), & 
							ez_g(k,j)/lmdd*t_int(1), er_g(k,j)/lmdd*t_int(1), & 
							vz_g(2,k,j), vr_g(2,k,j), vz_g(1,k,j), vr_g(1,k,j), &
							vtheta_g(1,k,j), jetheta_g(k,j), rho_g_SEE(1,k,j), jetheta_g_SEE(k,j)
				END DO
			END DO
		close(40)

		
		open(77,file='.\output\current.dat',position='append')
			write(77,"(1(1X,I7), 1(1X,F20.10))") iloop, Itheta, Itheta_SEE
		close(77)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! 离子入射到壁面上的能量、通量以及角度
!		open(111, file='.\output\boundary1.dat', status='replace')
!			write(111, "(A200)") 'VARIABLES = "I" "J" "Z" "R" "ION ENERGY" "ION FLUX" "ANGLE"'
!			j = gcb(1,3)
!			do i=gcb(1,1),gcb(1,2)
!				vi_abs = sqrt(vz_g(2,i,j)*vz_g(2,i,j)+vr_g(2,i,j)*vr_g(2,i,j))
!				if (vi_abs > 0.) then
!					WRITE(111,"(2(1X,I4), 3(1X,F15.10), 1(1X,E15.5), 1(1X,F15.10))") &
!						i, j, (cb(1,1)+i*dzr)*lmdd, cb(1,3)*lmdd, &
!						0.5*m(2)*vi_abs*vi_abs/E/mfactor(2), rho_g(2,i,j)*vi_abs, &
!						acos(vz_g(2,i,j)/vi_abs)
!				end if
!			end do
!		close(111)
		
!		open(222, file='.\output\boundary2.dat', status='replace')
!			write(222, "(A200)") 'VARIABLES = "I" "J" "Z" "R" "ION ENERGY" "ION FLUX" "ANGLE"'
!			i = gcb(4,1)
!			do j=gcb(4,4),gcb(4,3),-1
!				vi_abs = sqrt(vz_g(2,i,j)*vz_g(2,i,j)+vr_g(2,i,j)*vr_g(2,i,j))
!				if (vi_abs > 0.) then
!					WRITE(222,"(2(1X,I4), 3(1X,F15.10), 1(1X,E15.5), 1(1X,F15.10))") &
!						i, j, cb(4,1)*lmdd, (cb(4,3)+j*drr)*lmdd, &
!						0.5*m(2)*vi_abs*vi_abs/E/mfactor(2), rho_g(2,i,j)*vi_abs, &
!						acos(vr_g(2,i,j)/vi_abs)
!				end if
!			end do
!		close(222)
		
!		open(333, file='.\output\boundary3.dat', status='replace')
!			write(333, "(A200)") 'VARIABLES = "I" "J" "Z" "R" "ION ENERGY" "ION FLUX" "ANGLE"'
!			j = gcb(4,3)
!			do i=gcb(4,1),gcb(4,2)
!				vi_abs = sqrt(vz_g(2,i,j)*vz_g(2,i,j)+vr_g(2,i,j)*vr_g(2,i,j))
!				if (vi_abs > 0.) then
!					WRITE(333,"(2(1X,I4), 3(1X,F15.10), 1(1X,E15.5), 1(1X,F15.10))") &
!						i, j, (cb(4,1)+i*dzr)*lmdd, cb(4,3)*lmdd, &
!						0.5*m(2)*vi_abs*vi_abs/E/mfactor(2), rho_g(2,i,j)*vi_abs, &
!						acos(vz_g(2,i,j)/vi_abs)
!				end if
!			end do
!		close(333)
		
!		open(444, file='.\output\boundary4.dat', status='replace')
!			write(444, "(A200)") 'VARIABLES = "I" "J" "Z" "R" "ION ENERGY" "ION FLUX" "ANGLE"'
!			i = gcb(4,2)
!			do j=gcb(4,3),gcb(4,4)
!				vi_abs = sqrt(vz_g(2,i,j)*vz_g(2,i,j)+vr_g(2,i,j)*vr_g(2,i,j))
!				if (vi_abs > 0.) then
!					WRITE(444,"(2(1X,I4), 3(1X,F15.10), 1(1X,E15.5), 1(1X,F15.10))") &
!						i, j, cb(4,2)*lmdd, (cb(4,3)+j*drr)*lmdd, &
!						0.5*m(2)*vi_abs*vi_abs/E/mfactor(2), rho_g(2,i,j)*vi_abs, &
!						acos(vr_g(2,i,j)/vi_abs)
!				end if
!			end do
!		close(444)
		
!		open(555, file='.\output\boundary5.dat', status='replace')
!			write(555, "(A200)") 'VARIABLES = "I" "J" "Z" "R" "ION ENERGY" "ION FLUX" "ANGLE"'
!			j = gcb(3,3)
!			do i=gcb(3,1),gcb(3,2)
!				vi_abs = sqrt(vz_g(2,i,j)*vz_g(2,i,j)+vr_g(2,i,j)*vr_g(2,i,j))
!				if (vi_abs > 0.) then
!					WRITE(555,"(2(1X,I4), 3(1X,F15.10), 1(1X,E15.5), 1(1X,F15.10))") &
!						i, j, (cb(3,1)+i*dzr)*lmdd, cb(3,3)*lmdd, &
!						0.5*m(2)*vi_abs*vi_abs/E/mfactor(2), rho_g(2,i,j)*vi_abs, &
!						acos(vz_g(2,i,j)/vi_abs)
!				end if
!		close(555)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		rho_g1 = 0.
		rhovz_g1 = 0.
		rhovr_g1 = 0.
		rhovtheta_g1 = 0.
		rhoegy_g1 = 0.
		rhovt_g1 = 0.
		rhovt2_g1 = 0.
		ez_g1 = 0.
		er_g1 = 0.
		phi_g1 = 0.

		rhoegy_g1_normal = 0.
		rho_g1_SEE = 0.
		rhovtheta_g1_SEE = 0.

		IncNum = 0.
		EmitNum = 0.
	
	end if
	
	!对一些全局变量统一清零
	nout = 0
	
	if( mod(iloop,10)==0 ) then
		! 输出图形显示
!		neo = rho(1,:,:)
!		nio = rho(2,:,:)
!		call viewer(phi,0,nzt,0,nrt,"phi")
!		call viewer(neo,0,nzt,0,nrt,"ne")
!		call viewer(nio,0,nzt,0,nrt,"ni")

		open(1,file='.\output\number.dat',position='append')
			write(1,"(3(1X,I7),2(1X,F20.10))") iloop, num(1), num(2), minval(phi), ave_egy(1), ave_egy(2)
		close(1)
	end if

end subroutine output


!subroutine viewer(array,s1,b1,s2,b2,name)
!    use AVDef
!	implicit none
!	character(10) :: name
!	integer :: s1,b1,s2,b2
!	real :: array(s1:b1,s2:b2)
!	integer(4) status
!	call faglStartWatch(array, status)
!	call faglLBound(array, (0, 50), status)
!	call faglShow(array, status)
!	call faglName(array, name, status)
!	call faglUpdate(array, status)
!end subroutine viewer

! ++++++++++++++++++++++++++++++++++++++++++++++
!
! 2009-05-08���Ի��cmprg_gradient�Ĵ�������˸��ƣ�������һ������
! 
! 2009-05-15 ���£�
! ���Poisson����ʱ�����ݽڵ�λ�õĲ�ͬҪʹ�ò�ͬ����ɢ��ָ�ʽ���ر���ڱ߽�ڵ㣩�������Ҫ�Խڵ����ͺ�λ�������жϡ�
! ��������������������������жϣ�
!
! (1) ro_boundary: ��ʾ�ڵ�����.
! =1: �ڽڵ㣻=2: �߽�ڵ㣨�ǽǵ㣩��=4: �ڽǵ㣻=4/3: ��ǵ�
!
! (2) cmprg_gradient����ʾ�߽�ڵ㣨�ǽǵ㣩�����ĸ��߽磨���ҡ��ϡ��£����Լ��ǵ������������߽繹��.
! =2����߽磻=4���ұ߽磻=5���±߽磻=6���ϱ߽硣
! ע��1) ���ڽǵ��������߽�Ĺ����㣬���cmprg_gradient����ֵҪ���ӣ�
!     2) interior�߽��ϵĽڵ��൱���ڵ㣬��cmprg_gradient������1��
!
! (3) cmprg_bdtype����ʾ�߽�����
! cmprg_bdtypeֵ��ibֵ��ͬ���ڽǵ㴦��Ҫ���ӣ�������ӷ�����������ز��֣�
! �����ڽǵ��interior�߽磬cmprg_bdtype = 3��
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

	real lcb(num_a, 4), volumet(num_a), bt(0:IMAX1,0:IMAX2) ! lcb��ÿ���ߵĳ���
	real fluxrho_i, flux
	real krypton(energy_max,4) ! �������ײ���棬��2��3��4�зֱ�����ԡ������͵�����ײ����
	real, allocatable :: B(:,:)
	integer i, j, ii, jj, inum_a, k, insp, ipre, jpre, num_ea(num_a)
	integer inject_num, fnum_inject

	dzr = (z2-z1) / nzt ! ����ĺ���ߴ�
	drr = (r2-r1) / nrt ! ����ľ���ߴ�
	do i = 1,num_a
		dz(i) = (cb(i,2)-cb(i,1)) / nz(i) ! ��i���Ӽ������������ĺ���ߴ�
		dr(i) = (cb(i,4)-cb(i,3)) / nr(i) ! ��i���Ӽ������������ľ���ߴ�
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
	
	! ro_boundary���������interior����֮��ı߽紦���ܶ�ʱ��Ҫ��ϵ����ͬʱ��Ϊ����Possion���̵��ж�����֮һ
	! һ��أ�ro_boundary = 2;
	! �ر�أ����ڱ߽罻�㣬ro_boundary = 4 or 4/3.
	ro_boundary = 0.
	do inum_a = 1,num_a
	    ! ���ÿ���Ӽ���������������Ͻڵ㣨�������߽罻�㣩����ro_boundary
	    do i = 1,2
		    jpre = gcb(inum_a,i+2)      ! jpre����������inum_a���ϡ��±߽��������
		    if(ib(inum_a,i+2) /= 3) then        ! ��interior�߽�
		        do ipre = gcb(inum_a,1)+1,gcb(inum_a,2)-1       ! ipre��jpre�߽��ϵĽڵ㣨�������߽�˵㣩
		            ro_boundary(ipre,jpre) = 2.
		        end do
		    end if
		    ipre = gcb(inum_a,i)        ! ipre����������inum_a�����ұ߽��������
		    if(ib(inum_a,i) /= 3) then      ! ��interior�߽�
		        do jpre = gcb(inum_a,3)+1,gcb(inum_a,4)-1       ! jpre��ipre�߽��ϵĽڵ㣨�������߽�˵㣩
		            ro_boundary(ipre,jpre) = 2.
		        end do
		    end if
		end do
		
		! ���ÿ���Ӽ���������ĸ��߽罻�㴦��ro_boundary
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
    
    ! cmprg_gradient������߽�糡ʱ���ж�������ͬʱ��Ϊ����Poisson���̵��ж�����֮һ��������ʵ�ʼ��㡣
    ! ���������Ӽ�������ķ�interior���ͱ߽磬cmprg_gradient = 2, 4, 5, 6 for left, right, lower, upper boundaries, respectively.
    ! ����interior���ͱ߽�ڵ���ڽڵ㣬cmprg_gradient = 1.
    ! ע����һ��interior���ͱ߽�ǵ����ͬʱ����2���߽磬���cmprg_gradientҪ���ӡ�
    !     ���磬�ǵ㼴������߽磬���������±߽磬��cmprg_gradient = 2+5 = 7��
	cmprg_gradient = 0
	do i = 1,num_a
		do j = 1,4
			if(ib(i,j) /= 3) then       ! ��interior�߽�
				if(j < 3) then        ! �Ӽ�����������ұ߽�
					do k = 0,nr(i)
						ipre = gcb(i,j)
						jpre = gcb(i,3)+k
						if (ro_boundary(ipre,jpre) == 2) then        ! �Ǳ߽罻�������ڵ�
							if(j == 1) cmprg_gradient(ipre,jpre) = 2
							if(j == 2) cmprg_gradient(ipre,jpre) = 4
						else if (ro_boundary(ipre,jpre) /= 1) then      ! �߽罻��
							if(j == 1) cmprg_gradient(ipre,jpre) = cmprg_gradient(ipre,jpre)+2
							if(j == 2) cmprg_gradient(ipre,jpre) = cmprg_gradient(ipre,jpre)+4
						end if
					end do
				else        ! �Ӽ���������ϡ��±߽�
					do k = 0,nz(i)
						ipre = gcb(i,1)+k
						jpre = gcb(i,j)
						if(ro_boundary(ipre,jpre) == 2) then        ! �Ǳ߽罻�������ڵ�
							if(j == 3) cmprg_gradient(ipre,jpre) = 5
							if(j == 4) cmprg_gradient(ipre,jpre) = 6
						else if (ro_boundary(ipre,jpre) /= 1) then      ! �߽罻��
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
	    do j = 0,nz(i)        ! �Ӽ�������ĺ����ڽڵ�
	        do k = 0,nr(i)        ! �Ӽ�������ľ����ڽڵ�
				ipre = gcb(i,1)+j
				jpre = gcb(i,3)+k
				if (cmprg_gradient(ipre,jpre) == 0) cmprg_gradient(ipre,jpre) = 1
			end do
		end do
    end do
    ! +++++++++++END++++++++++++++++++++++++++++
    
    ! +++++++++++2009-05-15+++++++++++++++++++++
    ! cmprg_bdtype������Poisson���̵��ж�����֮һ����ʶ�˽ڵ����ڵı߽����ͣ�������ʵ�ʼ��㡣
    ! ���������Ӽ�������ķ�interior���ͱ߽磬cmprg_bdtype = ib(i,j)��
    ! ����interior���ͱ߽���ڽڵ㣬cmprg_bdtype = 3��
    ! ע����һ��interior���ͱ߽�Ķ˵����ͬʱ����2���ӱ߽磬���cmprg_bdtypeҪ���ӡ�
    ! �������ӣ����ڵ����ڵ�һ������dielectric surface�߽磬��һ������periodic�߽磬��cmprg_bdtype = 24 or 42����"2"��"4"����ϡ�
	cmprg_bdtype = 0
	do i = 1,num_a
		do j = 1,4
			if(ib(i,j) /= 3) then       ! ��interior�߽�
			    if(ib(i,j)==2 .OR. ib(i,j)==4 .OR. ib(i,j)==5) then     ! Ŀǰֻ�����˳���interior�߽�֮���2��4��5��߽�����
			        if(j < 3) then        ! �Ӽ�����������ұ߽�
					    do k = 0,nr(i)
						    ipre = gcb(i,j)
						    jpre = gcb(i,3)+k
						    if (ro_boundary(ipre,jpre) == 2) then        ! �Ǳ߽罻�������ڵ�
						        cmprg_bdtype(ipre,jpre) = ib(i,j)
						    else if (ro_boundary(ipre,jpre) /= 1) then      ! �߽罻��
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
				    else        ! �Ӽ���������ϡ��±߽�
					    do k = 0,nz(i)
						    ipre = gcb(i,1)+k
						    jpre = gcb(i,j)
						    if(ro_boundary(ipre,jpre) == 2) then        ! �Ǳ߽罻�������ڵ�
							    cmprg_bdtype(ipre,jpre) = ib(i,j)
						    else if (ro_boundary(ipre,jpre) /= 1) then      ! �߽罻��
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
	    do j = 0,nz(i)        ! �Ӽ�������ĺ����ڽڵ�
	        do k = 0,nr(i)        ! �Ӽ�������ľ����ڽڵ�
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
                ro_boundary_flag(i,j) = 1	! �ڵ�
            else if ( abs(ro_boundary(i,j)-2.0)<=float_error ) then
                ro_boundary_flag(i,j) = 2	! �߽��
            else if ( abs(ro_boundary(i,j)-4./3.)<=float_error ) then
                ro_boundary_flag(i,j) = 3	! ��ǵ�
            else if ( abs(ro_boundary(i,j)-4.0)<=float_error ) then
                ro_boundary_flag(i,j) = 4	! �ڽǵ�
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
	
	! ȷ�����ɷ��̱��ı�ı������°ݳ��Ⱥ͵�������Ƶ����Ӧ�ر��ı�
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
	octp = omigc(1)/omigp ! ���ӻ���Ƶ�ʺ͵�����Ƶ�ʵıȣ����ڵ����˶��ļ���

	t_parameter = 0.5*m*veth**2/e

	! parameter of ion
	vir_ion0 = SQRT(kb*t_int(1)*11600/m(2)) ! bohm velocity
	fluxqi = e*nero*vir_ion0*dt/OMIGP ! һ��ʱ�䲽�������ӵĵ�ɸ���
	
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
		! ����������ϵ
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

	! ����������ϵ
	volumet(:) = (cb(:,4)-cb(:,3)) *(cb(:,2)-cb(:,1)) *lmdd**2
!	do insp = 1, nsp
!		num(insp) = 5.e6 ! ����ϵͳ�ȶ�ʱÿ�����ӵ���ʵ����
!		! ����ÿ�������Ӵ�����ٸ���ʵ����
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
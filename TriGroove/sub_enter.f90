subroutine quasienter
	use global
	use discharge
	use grid
	use wall
	implicit none

	integer inum_a,insp,i,ipre,jpre,j,mc,num_in,n_cell_in,ipret,jpret
	integer ic1(2,imax1*imax2) ! 各个网格内的电子和离子数
	integer :: num_at_inlet(nsp,num_a,imax1,imax2)=0 ! 在入口每个网格节点内的电子和离子数
	integer :: net_at_boundary = 0
	integer :: net_at_cell(num_a,imax1,imax2) = 0
	real :: pro_at_cell(num_a,imax1,imax2) = 0.
	real fluxrho_i, belta, ranum, flux, a, v1, v2, v3, vhall, vir_ion0
	
	vhall = 2.e4/0.02
!	vhall = 0.
!	vir_ion0 = SQRT(kb*t_int(1)*11600/m(2))
	vir_ion0 = 0.


	insp = 2
	fluxrho_i = nero * sqrt(Kb*11600*t_int(insp)/(2*pi*m(insp))) ! injected net ion flux density

	do inum_a = 1, num_a
		do i = 1, 4
			if(ib(inum_a,i) == 5) then
				if(i==4) then
					flux = fluxrho_i * dt/omigp  / fnum(inum_a) &
								* (cb(inum_a,2)-cb(inum_a,1))*lmdd
					a = flux + num_inr(insp)
					num_in = int(a)
					num_inr(insp) = a - num_in
					do j = num(insp)+1, num(insp)+num_in
						call loadv(vz(insp,j), t_int(insp), insp)
						vz(insp,j) = vz(insp,j)/veth
100						call loadv(vr(insp,j), t_int(insp), insp)
						if (vr(insp,j) >= 0) then
							goto 100
						end if
						vr(insp,j) = vr(insp,j)/veth
!						vr(insp,j) = vr(insp,j) +(7-i*2)*vir_ion0/veth
						T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2)
!						if (insp == 1) then
!							call loadv(vtheta(j), t_int(insp), insp)
!							vtheta(j) = vtheta(j)/veth
!							vz(insp,j) = vz(insp,j) +vhall/veth
!							T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2+vtheta(j)**2)
!						end if
						r(insp,j) = cb(inum_a,i) +(7-i*2)*DR(INUM_A)*0.0001
						r_ini(insp,j) = r(insp,j)
						call random(ranum)
						z(insp,j) = cb(inum_a,1) +(cb(inum_a,2)-cb(inum_a,1))*ranum
					end do
					num(insp) = num(insp) + num_in
				else
					write(*,*) "error in file SUB_ENTER"
					stop
				end if
			end if
		end do
	end do
	
	n_cell_in = 2 ! 对临近边界的n_cell_in排网格采用准中性边界条件
	belta = sqrt(2*Kb*11600*t_int(1)/m(1))/veth ! 电子无量纲化后的热速度
	ic = 0
    do insp = 1, 2
        call INDEXM(insp)
        ic1(insp,:) = IC(2,:)
    end do

	do inum_a = 1,num_a
	    do i = 1, 4
	        if(ib(inum_a,i)==5) then
				if(i==4) then
					do jpre = nr(inum_a) - n_cell_in, nr(inum_a)
						net_at_cell = 0
						net_at_boundary = 0
						pro_at_cell = 0
						do ipre = 1, nz(inum_a)
							ipret = gcb(inum_a,1) + ipre
							jpret = gcb(inum_a,3) + jpre
							mc = (ipret-1)*imax2 + jpret
							! 统计入口边界附近网格内离子比电子多的数目net_at_cell
							net_at_cell(inum_a,ipre,jpre) = ic1(2,mc) - ic1(1,mc)
							net_at_boundary = net_at_boundary + net_at_cell(inum_a,ipre,jpre)
						end do
						if (net_at_boundary > 0) then ! >0，喷入电子，反之则不喷
							net_at_boundary = 0
							do ipre = 1, nz(inum_a)
								if (net_at_cell(inum_a,ipre,jpre) > 0) then
									pro_at_cell(inum_a,ipre,jpre) = real(net_at_cell(inum_a,ipre,jpre))
									net_at_boundary = net_at_boundary + net_at_cell(inum_a,ipre,jpre)
								else
									pro_at_cell(inum_a,ipre,jpre) = 0.
								end if
							end do
							pro_at_cell(inum_a,:,jpre) = pro_at_cell(inum_a,:,jpre)/real(net_at_boundary)
							do ipre = 2, nz(inum_a)
								pro_at_cell(inum_a,ipre,jpre) = pro_at_cell(inum_a,ipre,jpre) &
															+ pro_at_cell(inum_a,ipre-1,jpre)
							end do
							do j = num(1) + 1, num(1) + net_at_boundary
								call random(ranum)
								ipre = 1
								do while (ranum > pro_at_cell(inum_a, ipre, jpre) &
											.AND. ipre < nz(inum_a))
									ipre = ipre + 1
								end do

								call random(ranum)
								r(1,j) = cb(inum_a,3) +(jpre-ranum)*dr(inum_a)
!								r(1,j) = cb(inum_a,i) +(7-i*2)*DR(INUM_A)*0.0001
	                			call random(ranum)
	                			z(1,j) = cb(inum_a,1) +(ipre-ranum)*dz(inum_a)

								call loadv(vz(1,j), t_int(1), 1)
								vz(1,j) = vz(1,j)/veth +vhall/veth
200								call loadv(vr(1,j), t_int(1), 1)
								if (vr(1,j)>=0) then
									goto 200
								end if
								vr(1,j) = vr(1,j)/veth
								call loadv(vtheta(j), t_int(1), 1)
								vtheta(j) = vtheta(j)/veth
								T(1,j) = t_parameter(1)*(VZ(1,j)**2+VR(1,j)**2+Vtheta(j)**2)
							end do
							num(1) = num(1) +net_at_boundary
						end if
					end do
				else
					write(*,*) "error in file SUB_ENTER"
					stop
				end if
			end if
		end do
	end do

end subroutine quasienter

subroutine reflux
	use global
	use discharge
	use grid
	use wall
	implicit none

	integer inum_a, insp, i, j, num_in
	real ranum, flux, a, vhall, vir_ion0, fluxrho_i
	
	! vhall = 2.e4/0.01
	vhall = 0.
	! vir_ion0 = SQRT(kb*t_int(1)*11600/m(2))
	vir_ion0 = 0.
	
	insp = 2
	fluxrho_i = nero * sqrt(Kb*11600*t_int(insp)/(2*pi*m(insp))) ! injected net ion flux density
	! The net fluxes injected from the source are assumed temporally constant and equal for both species.
	! The emitted flux means injected plus refluxed particles.
	do inum_a = 1, num_a
		do i = 1, 4
			if(ib(inum_a,i) == 5) then
				do insp = 1, 2
					if(i==4 .OR. i==3) then
						flux = fluxrho_i * dt/omigp  / fnum(inum_a) &
								* (cb(inum_a,2)-cb(inum_a,1))*lmdd &
								+ nout(insp, inum_a, i)
!						flux = fluxrho_i * dt/omigp  / fnum(inum_a) &
!								* (cb(inum_a,2)-cb(inum_a,1))*lmdd
						a = flux + num_inr(insp)
						num_in = int(a)
						num_inr(insp) = a - num_in
						do j = num(insp)+1, num(insp)+num_in
							call loadv(vz(insp,j), t_int(insp), insp)
							vz(insp,j) = vz(insp,j)/veth
100							call loadv(vr(insp,j), t_int(insp), insp)
							if (i==4 .and. vr(insp,j)>=0) then
								goto 100
							else if (i==3 .and. vr(insp,j)<=0) then
								goto 100
							end if
							vr(insp,j) = vr(insp,j)/veth
							if (insp == 1) then
								call loadv(vtheta(j), t_int(insp), insp)
								vtheta(j) = vtheta(j)/veth
							end if
							if (insp == 1) then
								vz(insp,j) = vz(insp,j) +vhall/veth
								T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2+Vtheta(insp)**2)
							end if
							if (insp == 2) then
								vr(insp,j) = vr(insp,j) +(7-i*2)*vir_ion0/veth
								T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2)
							end if
							r(insp,j) = cb(inum_a,i) +(7-i*2)*DR(INUM_A)*0.0001
							call random(ranum)
							z(insp,j) = cb(inum_a,1) +(cb(inum_a,2)-cb(inum_a,1))*ranum
						end do
						num(insp) = num(insp) + num_in
					end if
				end do
			end if
		end do
	end do

end subroutine reflux
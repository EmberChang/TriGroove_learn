PROGRAM MAIN
	use global
	use neutral
	use wall
	use grid_para
	use mod_dadi
	USE DFPORT

	implicit none

	integer sloop,i,j,nql,iloop0,inspa
	real ave_egy,egy_total
	real elapsed_time
	character loop
	character(8) start_time, end_time

	start_time = CLOCK()
	
	sloop = 2.3e5
	state = 1
	iloop0 = 0

	write(*,*) 'input 0,1 for continuing, new calculation:'
	read (*,*) nql
	write(*,*) 'input 1,2 for state, 1 means simulate atoms, 2 means simulate plasma'
	read (*,*) state

	call triangle_geometry
	call initial

	if(nql == 0) then
		write(*,*) 'Read the restart file ......'
        if(state==1) then
            OPEN(4,FILE='.\output\case.res',FORM='UNFORMATTED')
                inspa = 3
                read(4) iloop0,num(inspa),z(inspa,1:num(inspa)),r(inspa,1:num(inspa)), &
                       r_ini(inspa,1:num(inspa)),vz(inspa,1:num(inspa)),vr(inspa,1:num(inspa)),ifield1
            CLOSE(4)
        else if(state==2) then
            OPEN (44,FILE='.\output\case_11.res',STATUS='OLD',FORM='UNFORMATTED')
                 READ(44) iloop0, num, ((z(i,j),j=1,num(i)),i=1,nsp), ((r(i,j),j=1,num(i)),i=1,nsp), & 
                       ((r_ini(i,j),j=1,num(i)),i=1,nsp), ((vz(i,j),j=1,num(i)),i=1,nsp), & 
                       ((vr(i,j),j=1,num(i)),i=1,nsp), vtheta(1:num(1)), ((t(i,j),j=1,num(i)),i=1,nsp), & 
                       fnuma(1:num(nsp)), fnum, collectq, rho, conductor_charge, ifield2, pfactor, &
                       dt, lmdd, omigp, nout, ez, er
            CLOSE(44)
			ifield2 = 0
        end if
	else
		open(1,file='.\output\number.dat')
    		write(1,"(A160)") ' VARIABLES = "step" "num_e" "num_i" "phiwall" "egye" "egyi" '
		close(1)

		open(11,file='.\output\current.dat')
    		write(11,"(A60)") ' VARIABLES = "step" "Iez" "Iez_SEE" '
		close(11)

		open(110,file='.\output\current2.dat')
    		write(11,"(A60)") ' VARIABLES = "step" "Iez" '
		close(110)

		if(performance_diag==1) then
			open(111,file='.\output\performance.dat')
    			write(111,"(A260)") ' VARIABLES = "step" "time" "num_a" "Thrust" "Isp" "Ii" &
							"anode_current" "efficiency_discharge" "efficiency_utilization" &
							"efficiency_acceleration" "num_e" "num_i" "num_a1" "maxphi" &
							"conductor_charge" "phi_outter" "maxroe" "maxroi" "maxs" "nero" "ionrate" '
			close(111)
		end if

		if(iloop==1111) then
			open(1111,file='.\output\neutral.dat')
				write(1,"(A160)") 'variables = "step" "time" "nflux"'
			close(1111)
		end if
	end if

	elapsed_time = TIMEF()

100	do ILOOP = iloop0+1, sloop
		if(iloop>iloopa .AND. state==1) state = 2
		if(state==1) then

		else
			call move ! move必须放在dadi前，不然会出现不稳定
			call quasienter
!			call reflux
			call gather
!			do i = 1, 2
!				call smooth(rho(1,:,:), imax1, imax2)
!				call smooth(rho(2,:,:), imax1, imax2)
!			end do
			call dadi
			call output
			
			egy_total = 0.
			do j = 1,num(2)
				egy_total = egy_total + t(2,j)
			end do
			ave_egy = 0
			if(num(2)>0) ave_egy = egy_total/num(2)

			if(mod(iloop,10)==0 .and. nsp==1) then
				WRITE(*,"('iloop=', I6, 1X, 'electron=', I7, 1X, 'ave_egy=', F10.4)") &
									iloop,				num(1),				ave_egy
			else if (mod(iloop,10)==0 .and. nsp==2) then
				WRITE(*,"('iloop=', I6, 1X, 'electron=', I7, 1X, 'ion=', I7, 1X,' ave_egy=', F10.4)") &
									iloop,				num(1),			num(2),				ave_egy
			else if (mod(iloop,10)==0 .and. nsp==3) then
			    inspa=3
				WRITE(*,"('iloop=', I6, 1X, 'electron=', I7, 1X, 'ion=', I7, 1X,' atom=', I7, 1X, 'ave_egy=', F10.4)")  &
									iloop,				num(1),			num(2),			num(inspa),			ave_egy
			end if
		end if
	END DO

	end_time = CLOCK()
	write (*,*) 'The begin time is ', start_time
	write (*,*) 'The finish time is ', end_time

	write(*,*) "whether continue ? (y/n)"
	read(*,*) loop
	if (loop == 'y') then
		iloop0 = sloop
		write(*,*) "please input the further loop number"
		read(*,*) iloop
		sloop = sloop + iloop
		goto 100
	end if
	
END PROGRAM MAIN
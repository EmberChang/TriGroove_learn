subroutine diagnostic
	use constant
	use global
	use mcc_para
	use grid
	use wall
	use diag
	implicit none
	
	real zepre,repre,thetapre ! location of nonominized value
	real para(nsp,NumPar_MAX)
	integer i, j
	integer :: insp = 1
	character a, b
	integer ifile
	
	do ifile = 1, ipc_num
		if(ipc(ifile) == 0) then
			zepre = z(insp,ifile)
			repre = r(insp,ifile) + r0
			thetapre = 0
			if(ifile.le.9) then
				a = char(ifile+48)
				open(1,file='.\output\trace'//a//'.dat',position='append')
			else
				a = char(int(ifile/10)+48)
				b = char(mod(ifile,10)+48)
				open(1,file='.\output\trace'//a//b//'.dat', position='append')
			end if
			
				write(1,"(3(1X,F12.6))") zepre*lmdd, repre*lmdd, thetapre !LOCATION IN Z-R-PHI COORDINATE
			close(1)
		end if
	end do

end subroutine diagnostic


subroutine ave_free_path(insp,j)
	use diag
	use global
	use magnetic
	implicit none

	integer insp,j,init_ipre,init_jpre,last_ipre,last_jpre
	real init_period,last_period,time_free,init_cycles,last_cycles

	init_ipre = INT((init_z(j) - z1) / dzr)
	init_jpre = INT((init_r(j) - r1) / drr)
	last_ipre = INT((z(insp,j) - z1) / dzr)
	last_jpre = INT((r(insp,j) - r1) / drr)
	init_period = 1/(e*sqrt(bzt(init_ipre,init_jpre)**2+brt(init_ipre,init_jpre)**2)/(m(insp)*twopi))
	last_period = 1/(e*sqrt(bzt(last_ipre,last_jpre)**2+brt(last_ipre,last_jpre)**2)/(m(insp)*twopi))
	time_free = collision_time(j)/omigp*dt
	init_cycles = time_free/init_period
	last_cycles = time_free/last_period

	open(1,file='.\output\free_path.dat',position='append')
		write(1,"(4(1X,F10.5),(1X,E15.5),2(1X,F10.5),(1X,I5))") init_z(j)*lmdd,init_r(j)*lmdd,z(insp,j)*lmdd, &
				r(insp,j)*lmdd,time_free,init_cycles,last_cycles,collision_time(j)
	close(1)

	collision_time(j) = 0
	init_z(j) = z(insp,j)
	init_r(j) = r(insp,j)
end subroutine ave_free_path

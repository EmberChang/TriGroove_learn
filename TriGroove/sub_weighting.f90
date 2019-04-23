subroutine co_gtop(insp,j,ipret,jpret,a1,a2,a3,a4)
	use global
	implicit none
	integer insp, j, ipret, jpret
	real a1, a2, a3, a4, ztem, rtem
	real z_i, r_j

	ztem = -1.
	rtem = -1.

	do while (ztem<0. .or. ztem>1.)
		z_i		 =	z1 + ipret * dzr
		ztem	 =	(z(insp,j) - z_i) / dzr
		if (ztem>1.) then
			ipret = ipret +1
			IF(ipret == nzt) then
				ipret = ipret-1
				z(insp,j) = z2 - dzr*0.001
			end if
		else if (ztem<0.) then
			ipret = ipret -1
			IF(ipret == -1) then
				ipret = ipret+1
				z(insp,j) = z1 + dzr*0.001
			end if
		end if
	end do

	do while (rtem<0. .or. rtem>1.)
		r_j		 =	r1 + jpret * drr
		rtem	 =	(r(insp,j) - r_j) / drr
		if (rtem>1.) then
			jpret = jpret +1
			IF(jpret == nrt) then
				jpret = jpret-1
				r(insp,j) = r2 - drr*0.001
			end if
		else if (rtem<0.) then
			jpret = jpret -1
			IF(jpret == -1) then
				jpret = jpret+1
				r(insp,j) = r1 + drr*0.001
			end if
		end if
	end do

	a1 = (1. - rtem) * (1. - ztem)					!i,		j
	a2 = rtem * (1. - ztem)							!i,		j+1
	a3 = (1. - rtem) * ztem							!i+1,	j
	a4 = rtem * ztem								!i+1,	j+1
end subroutine co_gtop


subroutine gtop(insp,j,parapre,para,imax10,imax20,ipret,jpret,a1,a2,a3,a4) !grid to particle
	use global
	use magnetic
	implicit none
	integer :: imax10, imax20
	real :: para(0:imax10, 0:imax20)
	real a1, a2, a3, a4, parapre
	integer i, j, ipret, jpret, inum_a, insp

	if(a2 * a4 == 0) then
		parapre = para(ipret, jpret) * a1 + para(ipret+1, jpret) * a3
	else if (a3 * a4 == 0) then
		parapre = para(ipret, jpret) * a1 + para(ipret, jpret+1) * a2
	else
		parapre = para(ipret, jpret) * a1 + para(ipret, jpret+1) * a2 &
				+ para(ipret+1, jpret) * a3 + para(ipret+1, jpret+1) * a4
	end if
end subroutine gtop


subroutine co_ptogt(insp,j,inum_a,ipret,jpret,p1,p2,p3,p4)
	use global
	use grid
	use neutral 
	use neutral
	implicit none
	integer insp,j,inum_a,ipret,jpret
	real ztem,rtem,r0tor,a1,a2,a3,a4,p1,p2,p3,p4
	real z_i, r_j, r_jplus1

	ztem = -1.
	rtem = -1.

	do while (ztem<0. .or. ztem>1.)
		z_i		 =	z1 + ipret * dzr
		ztem	 =	(z(insp,j) - z_i) / dzr
		if (ztem>1.) then
			ipret = ipret +1
		else if (ztem<0.) then
			ipret = ipret -1
		end if
	end do

	do while (rtem<0. .or. rtem>1.)
		r_j		 =	r1 + jpret * drr
		rtem	 =	(r(insp,j) - r_j) / drr
		if (rtem>1.) then
			jpret = jpret +1
		else if (rtem<0.) then
			jpret = jpret -1
		end if
	end do

!	z_i		 =	z1 + ipret * dzr
!	ztem	 =	(z(insp,j) - z_i) / dzr
!	r_j		 =	r1 + jpret * drr
!	rtem	 =	(r(insp,j) - r_j) / drr
		
	a1 = (1. - rtem) * (1. - ztem)					!i,		j
	a2 = (1. - rtem) * ztem							!i+1,	j
	a3 = rtem * (1. - ztem)							!i,		j+1
	a4 = rtem * ztem								!i+1,	j+1
	p1 = fnum(inum_a)*a1/volume_cell(jpret)
	p2 = fnum(inum_a)*a2/volume_cell(jpret)
	p3 = fnum(inum_a)*a3/volume_cell(jpret+1)
	p4 = fnum(inum_a)*a4/volume_cell(jpret+1)
	if(insp == 3) then
		p1 = p1*fnuma(j)
		p2 = p2*fnuma(j)
		p3 = p3*fnuma(j)
		p4 = p4*fnuma(j)
	end if	  
end subroutine co_ptogt


subroutine ptogt(density,nsp1,imax10,imax20,j,insp,para,ipret,jpret,p1,p2,p3,p4)
	use constant
	use global
	use grid
	use neutral
	use grid_para
	implicit none
	integer :: imax10,imax20,nsp1
	real :: density(nsp1,0:imax10,0:imax20)
	real p1,p2,p3,p4
	integer i,j,k,ipret,jpret,insp,inum_a
	real para
	density(insp, ipret, jpret) = density(insp, ipret, jpret) + p1*para*ro_boundary(ipret, jpret)
	density(insp, ipret+1, jpret) = density(insp, ipret+1, jpret) + p2*para*ro_boundary(ipret+1,jpret)
	density(insp, ipret, jpret+1) = density(insp, ipret, jpret+1) + p3*para*ro_boundary(ipret,jpret+1)
	density(insp, ipret+1, jpret+1) = density(insp, ipret+1, jpret+1) + p4*para*ro_boundary(ipret+1,jpret+1)
end subroutine ptogt


subroutine linear_weighting(insp,inum_a,ibound,k,j,para,num_a1,nsp1,n1,imax10,para1,flag)

	! 计算沉积到边界上的粒子数、电荷数或者能量
	use global
	use grid_para
	implicit none

	integer nsp1, num_a1, n1, imax10
	integer insp, inum_a, ibound, k, j, flag, jpre, ipre
	real para(nsp1, num_a1, n1, 0:imax10)
	real ztem, rtem, para1

	if(ibound == 1) then
		jpre = int((r(insp,j) - cb(inum_a,3)) / dr(inum_a))
		if (jpre >= nr(inum_a)) jpre = nr(inum_a) -1
		rtem = (r(insp,j)-cb(inum_a,3)-jpre*dr(inum_a)) / dr(inum_a)
		if(jpre == 0) then
			para(insp,inum_a,k,jpre) = para(insp,inum_a,k,jpre) + 2*(3-2*flag)*para1*fnum(inum_a)*(1-rtem)
			para(insp,inum_a,k,jpre+1) = para(insp,inum_a,k,jpre+1) + (3-2*flag)*para1*fnum(inum_a)*rtem
		else if(jpre == (nr(inum_a)-1)) then
			para(insp,inum_a,k,jpre) = para(insp,inum_a,k,jpre) + (3-2*flag)*para1*fnum(inum_a)*(1-rtem)
			para(insp,inum_a,k,jpre+1) = para(insp,inum_a,k,jpre+1) + 2*(3-2*flag)*para1*fnum(inum_a)*rtem
		else
			para(insp,inum_a,k,jpre) = para(insp,inum_a,k,jpre) + (3-2*flag)*para1*fnum(inum_a)*(1-rtem)
			para(insp,inum_a,k,jpre+1)=para(insp,inum_a,k,jpre+1) + (3-2*flag)*para1*fnum(inum_a)*rtem
		end if
	else
		ipre = int((z(insp,j) - cb(inum_a,1)) / dz(inum_a))
		if (ipre >= nz(inum_a)) ipre = nz(inum_a) -1
		ztem = (z(insp,j)-cb(inum_a,1)-ipre*dz(inum_a)) / dz(inum_a)
		if(ipre == 0) then
			para(insp,inum_a,k+2,ipre) = para(insp,inum_a,k+2,ipre) + 2*(3-2*flag)*para1*fnum(inum_a)*(1-ztem)
			para(insp,inum_a,k+2,ipre+1) = para(insp,inum_a,k+2,ipre+1) + (3-2*flag)*para1*fnum(inum_a)*ztem
		else if(ipre == (nz(inum_a)-1)) then
			para(insp,inum_a,k+2,ipre) = para(insp,inum_a,k+2,ipre) + (3-2*flag)*para1*fnum(inum_a)*(1-ztem)
			para(insp,inum_a,k+2,ipre+1) = para(insp,inum_a,k+2,ipre+1) + 2*(3-2*flag)*para1*fnum(inum_a)*ztem
		else
			para(insp,inum_a,k+2,ipre) = para(insp,inum_a,k+2,ipre) + (3-2*flag)*para1*fnum(inum_a)*(1-ztem)
			para(insp,inum_a,k+2,ipre+1) = para(insp,inum_a,k+2,ipre+1) + (3-2*flag)*para1*fnum(inum_a)*ztem
		end if
	end if
end subroutine linear_weighting
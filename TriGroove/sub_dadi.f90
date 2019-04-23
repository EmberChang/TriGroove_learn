module poisson
	use constant
	use global
	use ion
	use wall
	use mod_dadi
	use grid_para

	implicit none

	integer, parameter :: DADI_DEBUG = 0
	integer :: init_field_flag = 1
	real, parameter :: tol_test = 1.e-3
	real :: dirichletb(0:imax1,0:imax2) = 1.
	real :: rho_normalized(nsp,0:imax1,0:imax2) = 0.
	real u(0:imax1,0:imax2), uwork(0:imax1,0:imax2), ustor(0:imax1,0:imax2), ustar(0:imax1,0:imax2)
	real a_xgeom(0:imax1,0:imax2), b_xgeom(0:imax1,0:imax2), c_xgeom(0:imax1,0:imax2)
	real a_ygeom(0:imax1,0:imax2), b_ygeom(0:imax1,0:imax2), c_ygeom(0:imax1,0:imax2)
	real source(0:imax1,0:imax2), serho(0:imax1,0:imax2)
	real u_in(0:imax1,0:imax2)
	real dx, dy
	real :: dadi_dt = 0.0
end module poisson


subroutine dadi
	use poisson
	use discharge
	implicit none

	integer i,j,k,ipre,jpre,ncx,ncy
	integer itermax
	
	ncx = imax1
	ncy = imax2
	dx = dzr
	dy = drr
	
	itermax = 100
	
	call get_charge_density(ncx, ncy)
    
    if(init_field_flag) then
		dadi_dt = 0.1*(dx*dx+dy*dy)
		call init_dadi_arrays(nzt,nrt) ! 初始化泊松方程求解器
	    u = 0
	    u_in = 0
	    init_field_flag = 0
    else
        u = u_in
    end if

    call dadisolve(itermax,ncx,ncy)
	
	phi = u

	! neumann
!	do i = 0, ncx
!		phi(i,:) = phi(i,:) - phi(i,ncy)
!	end do

    do i = 0, ncx
        do j = 0, ncy
            if(cmprg_bdtype(i,j) == 0) phi(i,j) = 0.
        end do
    end do

    call gradient(ncx, ncy)
  
end subroutine dadi


subroutine get_charge_density(ncx,ncy)
	use poisson
	use grid_para
	use ion
	implicit none
	
	real sss,vir_ion
	integer i,j,k,ncx,ncy,ipre,jpre,insp,inum_a1,inum_a2
	real conductor_phi
	
	source = 0.0
	rho_normalized = rho /nero
	conductor_phi = conductor_charge /(capacitance *t_int(1))
	
	if (nsp /= 2) then
		write (*, "('Error: Wrong Particle Species Number: ', I4)") nsp
		stop
	end if

	do insp = 1, nsp
		source(0:ncx, 0:ncy) = source(0:ncx, 0:ncy) + rho_normalized(insp,0:ncx,0:ncy)*q1(insp)
	end do

  	inum_a1 = 1
	inum_a2 = num_a
	do i = inum_a1, inum_a2
		do j = 1, 4
		    if(ib(i,j) == 2) then
				if(j < 3) then
					do k = 0, nr(i)
						ipre = gcb(i,j)
						jpre = gcb(i,3) +k
						serho(ipre,jpre) = collectq(i,j,k)

						if (ro_boundary_flag(ipre,jpre)==2) then
							source(ipre,jpre) = source(ipre,jpre) + collectq(i,j,k)/dx
						else if (ro_boundary_flag(ipre,jpre)==3) then
							if (cmprg_bdtype(ipre,jpre)==22) then
								source(ipre,jpre) = source(ipre,jpre) + collectq(i,j,k)/dx/3
							else
								write (*, "('Error: Exception on boundary value ', I8)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						else if (ro_boundary_flag(ipre,jpre)==4) then
							if (cmprg_bdtype(ipre,jpre)==22 .OR. cmprg_bdtype(ipre,jpre)==42 &
								.OR. cmprg_bdtype(ipre,jpre)==24) then
								source(ipre,jpre) = source(ipre,jpre) + collectq(i,j,k)/dx
							else
								write (*, "('Error: Exception on boundary value ', I8)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						else
							write (*, "('Error: Exception on RO_BOUNDARY value ', I8)") ro_boundary_flag(ipre,jpre)
							stop
						end if
						
					end do
				else
					do k = 0,nz(i)
						ipre = gcb(i,1) +k
						jpre = gcb(i,j)
						serho(ipre,jpre) = collectq(i,j,k)

						if (ro_boundary_flag(ipre,jpre)==2) then
							source(ipre,jpre) = source(ipre,jpre) + collectq(i,j,k)/dy
						else if (ro_boundary_flag(ipre,jpre)==3) then
							if (cmprg_bdtype(ipre,jpre)==22) then
								source(ipre,jpre) = source(ipre,jpre) + collectq(i,j,k)/dy/3
							else
								write (*, "('Error: Exception on boundary value ', I8)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						else if (ro_boundary_flag(ipre,jpre)==4) then
							if (cmprg_bdtype(ipre,jpre)==22 .OR. cmprg_bdtype(ipre,jpre)==42 &
								.OR. cmprg_bdtype(ipre,jpre)==24) then
								source(ipre,jpre) = source(ipre,jpre) + collectq(i,j,k)/dy
							else
								write (*, "('Error: Exception on boundary value ', I8)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						else
							write (*, "('Error: Exception on RO_BOUNDARY value ', I8)") ro_boundary_flag(ipre,jpre)
							stop
						end if
						
					end do
				end if

			else if(ib(i,j) == 5) then
				if(j < 3) then
					do k = 0, nr(i)
						ipre = gcb(i,j)
						jpre = gcb(i,3) +k
						dirichletb(ipre,jpre) = 0 /t_int(1)
					end do
				else
					do k = 0,nz(i)
						ipre = gcb(i,1) +k
						jpre = gcb(i,j)
						dirichletb(ipre,jpre) = 0 /t_int(1)
					end do
				end if

!			else if(ib(i,j) == 7) then
!				if(j < 3) then
!					do k = 0, nr(i)
!						ipre = gcb(i,j)
!						jpre = gcb(i,3) +k
!						dirichletb(ipre,jpre) = ud /t_int(1)
!					end do
!				else
!					do k = 0,nz(i)
!						ipre = gcb(i,1) +k
!						jpre = gcb(i,j)
!						dirichletb(ipre,jpre) = ud /t_int(1)
!					end do
!				end if

			else if(ib(i,j) == 11) then
				if(j < 3) then
					do k = 0, nr(i)
						ipre = gcb(i,j)
						jpre = gcb(i,3) +k
						dirichletb(ipre,jpre) = conductor_phi
					end do
				else
					do k = 0,nz(i)
						ipre = gcb(i,1) +k
						jpre = gcb(i,j)
						dirichletb(ipre,jpre) = conductor_phi
					end do
				end if

			end if

		end do
	end do

	! 周期性条件的人为强制
	source(0,:) = (source(0,:) + source(ncx,:))/2
	source(ncx,:) = source(0,:)

!	if( mod(iloop,10)==0 ) then
!		! 输出图形显示
!		call viewer(source, 0, imax1, 0, imax2, "source")
!		call viewer(serho, 0, imax1, 0, imax2, "serho")
!	end if
    
end subroutine get_charge_density


subroutine init_dadi_arrays(ncx,ncy)
	use poisson
	use grid_para
	implicit none
	
	integer i,j,k,ipre,jpre,ncx,ncy,insp,inum_a1,inum_a2
	integer rho_flag, gradient_flag, bdtype_flag

    a_xgeom = 0.
    c_xgeom = 0.
    b_xgeom = 0.
    a_ygeom = 0.
    c_ygeom = 0.
    b_ygeom = 0.

    ! Fixing coefficients
    do j = 0,ncy
        do i = 0,ncx
            rho_flag = ro_boundary_flag(i,j)
            gradient_flag = cmprg_gradient(i,j)
            bdtype_flag = cmprg_bdtype(i,j)

            select case (rho_flag)
            case (1)
                a_xgeom(i,j) = 1./dx/dx
                c_xgeom(i,j) = 1./dx/dx
                b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                a_ygeom(i,j) = 1./dy/dy
                c_ygeom(i,j) = 1./dy/dy
                b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)

            case (2)
				if (bdtype_flag == 2) then
                    if (gradient_flag == 2) then ! left
                        a_xgeom(i,j) = 0.
                        c_xgeom(i,j) = 2./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 1./dy/dy
                        c_ygeom(i,j) = 1./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 4) then ! right
                        a_xgeom(i,j) = 2./dx/dx
                        c_xgeom(i,j) = 0.
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 1./dy/dy
                        c_ygeom(i,j) = 1./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 5) then ! down
                        a_xgeom(i,j) = 1./dx/dx
                        c_xgeom(i,j) = 1./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 0.
                        c_ygeom(i,j) = 2./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 6) then ! up
                        a_xgeom(i,j) = 1./dx/dx
                        c_xgeom(i,j) = 1./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 2./dy/dy
                        c_ygeom(i,j) = 0.
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    end if
                else if (bdtype_flag == 4) then
					if (gradient_flag == 2) then ! left
						a_xgeom(i,j) = 1./dx/dx
						c_xgeom(i,j) = 1./dx/dx
						b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
						a_ygeom(i,j) = 1./dy/dy
						c_ygeom(i,j) = 1./dy/dy
						b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
					else if (gradient_flag == 4) then ! right
                        a_xgeom(i,j) = 1./dx/dx
                        c_xgeom(i,j) = 1./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 1./dy/dy
                        c_ygeom(i,j) = 1./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
					else
						write (*, "('Error: Exception of Periodic Boundary')")
						stop
					end if
				else if (bdtype_flag == 5) then
					if (gradient_flag == 6) then
						! neumann
						a_xgeom(i,j) = 0.
						c_xgeom(i,j) = 0.
						b_xgeom(i,j) = 1.
						a_ygeom(i,j) = 0.
						c_ygeom(i,j) = 0.
						b_ygeom(i,j) = 1.
					else
						write (*, "('Error: Exception on boundary type ', I8)") bdtype_flag
                        stop
					end if
				else
                    write (*, "('Error: Undefined Boundary Type ', I8)") bdtype_flag
					stop
                end if

            case (3)
                if (bdtype_flag == 22) then
                    if (gradient_flag == 7) then
						a_xgeom(i,j) = 2./3./dx/dx
                        c_xgeom(i,j) = 4./3./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 2./3./dy/dy
                        c_ygeom(i,j) = 4./3./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 9) then
						a_xgeom(i,j) = 4./3./dx/dx
                        c_xgeom(i,j) = 2./3./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 2./3./dy/dy
                        c_ygeom(i,j) = 4./3./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 10) then
						a_xgeom(i,j) = 4./3./dx/dx
                        c_xgeom(i,j) = 2./3./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 4./3./dy/dy
                        c_ygeom(i,j) = 2./3./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 8) then
						a_xgeom(i,j) = 2./3./dx/dx
                        c_xgeom(i,j) = 4./3./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 4./3./dy/dy
                        c_ygeom(i,j) = 2./3./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    end if
                else
                    write (*, "('Error: Exception on 外角点边界 with boundary type ', I8)") bdtype_flag
				    stop
                end if

            case (4)
                if (bdtype_flag == 22) then
                    if (gradient_flag == 7) then
						a_xgeom(i,j) = 0.
                        c_xgeom(i,j) = 2./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 0.
                        c_ygeom(i,j) = 2./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 9) then
						a_xgeom(i,j) = 2./dx/dx
                        c_xgeom(i,j) = 0.
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 0.
                        c_ygeom(i,j) = 2./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 10) then
						a_xgeom(i,j) = 2./dx/dx
                        c_xgeom(i,j) = 0.
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 2./dy/dy
                        c_ygeom(i,j) = 0.
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else if (gradient_flag == 8) then
						a_xgeom(i,j) = 0.
                        c_xgeom(i,j) = 2./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 2./dy/dy
                        c_ygeom(i,j) = 0.
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    end if
                else if (bdtype_flag == 24 .OR. bdtype_flag == 42) then
                    if (gradient_flag == 7 .OR. gradient_flag == 9) then
						a_xgeom(i,j) = 1./dx/dx
                        c_xgeom(i,j) = 1./dx/dx
                        b_xgeom(i,j) = a_xgeom(i,j) +c_xgeom(i,j)
                        a_ygeom(i,j) = 0.
                        c_ygeom(i,j) = 2./dy/dy
                        b_ygeom(i,j) = a_ygeom(i,j) +c_ygeom(i,j)
                    else
                        write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") bdtype_flag
                        stop
                    end if
                else if (bdtype_flag == 45 .OR. bdtype_flag == 54) then
                    if (gradient_flag == 8 .OR. gradient_flag == 10) then
						a_xgeom(i,j) = 0.
						c_xgeom(i,j) = 0.
						b_xgeom(i,j) = 1.
						a_ygeom(i,j) = 0.
						c_ygeom(i,j) = 0.
						b_ygeom(i,j) = 1.
                    else
                        write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") bdtype_flag
                        stop
                    end if
                else
                    write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") bdtype_flag
				    stop
                end if

            case (0)
                continue

            case default
                write (*, "('Error: Exception on RO_BOUNDARY value ', F8.3)") ro_boundary(i,j)
                stop
            end select
        end do
    end do

!	call viewer(a_xgeom, 0, imax1, 0, imax2, "a_xgeom")
!	call viewer(b_xgeom, 0, imax1, 0, imax2, "b_xgeom")
!	call viewer(c_xgeom, 0, imax1, 0, imax2, "c_xgeom")
!	call viewer(a_ygeom, 0, imax1, 0, imax2, "a_ygeom")
!	call viewer(b_ygeom, 0, imax1, 0, imax2, "b_ygeom")
!	call viewer(c_ygeom, 0, imax1, 0, imax2, "c_ygeom")

end subroutine init_dadi_arrays


subroutine dadisolve(itermax,ncx,ncy)
	use poisson
	use discharge
	use grid
	use neutral
	use grid_para
	implicit none
	
	integer itermax, ncx, ncy
	integer i, j, ip1, im1, jp1, jm1
	integer iter, ndiscard
	real del_t, del_td, tptop, tpbot, ratio
	real rnorm, rsum, res, errchk, dxdxutrm, dydyutrm
	real iter_ave, res_ave

if (DADI_DEBUG==1) then
  	iter_ave = 0.0
  	res_ave = 0.0
end if
  
	! Start initial step-size as previous step size divided by n.
	del_t = dadi_dt/16.0
	del_td = 2.0*del_t
	ndiscard = 0

	! Residual normalization.
	rnorm = 0.0
	do j = 0, ncy
		do i = 0, ncx
			if (cmprg_bdtype(i,j)>0 .AND. cmprg_bdtype(i,j)/=5 &
				.AND. cmprg_bdtype(i,j)/=45 .AND. cmprg_bdtype(i,j)/=54) then
!				im1 = i-1
!				ip1 = i+1
!				jm1 = j-1
!				jp1 = j+1
!				if (im1 < 0) im1 = 0
!				if (ip1 > ncx) ip1 = ncx
!				if (jm1 < 0) jm1 = 0
!				if (jp1 > ncy) jp1 = ncy
!				! Residual terms.
!				dxdxutrm = a_xgeom(i,j)*u_in(im1,j) -b_xgeom(i,j)*u_in(i,j) +c_xgeom(i,j)*u_in(ip1,j)
!				dydyutrm = a_ygeom(i,j)*u_in(i,jm1) -b_ygeom(i,j)*u_in(i,j) +c_ygeom(i,j)*u_in(i,jp1)
!				! Residual sums.  Only include points outside of structures.
!				errchk = dxdxutrm +dydyutrm +source(i,j)
!				rnorm = rnorm +errchk*errchk
				rnorm = rnorm +source(i,j)*source(i,j)
			end if
		end do
	end do
	
	if (ABS(rnorm) < 1e-30) rnorm = 1.0
	rnorm = sqrt(rnorm)
	
	! Begin iteration.
	iter = 0
	res = 1.
	do while (iter<=itermax .AND. res > tol_test)
		iter = iter+1
		
		! Copy u into the work array and storage array.
		uwork = u
		ustor = u
		
		! Two advances of u via ADI at del_t.
		call adi(u, del_t, ncx, ncy)
		call adi(u, del_t, ncx, ncy)
		
		! One advance of uwork via ADI at 2*del_t.
		call adi(uwork, del_td, ncx, ncy)

!		call viewer(u,0,imax1,0,imax2,"u")
!		call viewer(uwork,0,imax1,0,imax2,"uwork")
!		call viewer(u-uwork,0,imax1,0,imax2,"u-uwork")
!		call viewer(u-ustor,0,imax1,0,imax2,"u-ustor")
	   
		! Calculate test parameter and normalized error.
		! For Dirichlet BCs, no need to worry about boundary
		! points since u, uwork, and ustor should be the same.
 
		tptop = 0.0
		tpbot = 0.0
		rsum = 0.0
		do j = 0, ncy
			do i = 0, ncx
				if (cmprg_bdtype(i,j)>0 .AND. cmprg_bdtype(i,j)/=5 &
					.AND. cmprg_bdtype(i,j)/=45 .AND. cmprg_bdtype(i,j)/=54) then
					! Test paramter sums.
					tptop = tptop + (u(i,j)-uwork(i,j))*(u(i,j)-uwork(i,j))
					tpbot = tpbot + (u(i,j)-ustor(i,j))*(u(i,j)-ustor(i,j))
					
					im1 = i-1
					ip1 = i+1
					jm1 = j-1
					jp1 = j+1
					if (im1 < 0) im1 = ncx-1
					if (ip1 > ncx) ip1 = 1
					if (jm1 < 0) jm1 = 0
					if (jp1 > ncy) jp1 = ncy
					! Residual terms.
					dxdxutrm = a_xgeom(i,j)*u(im1,j) -b_xgeom(i,j)*u(i,j) +c_xgeom(i,j)*u(ip1,j)
					dydyutrm = a_ygeom(i,j)*u(i,jm1) -b_ygeom(i,j)*u(i,j) +c_ygeom(i,j)*u(i,jp1)
					! Residual sums.  Only include points outside of structures.
					errchk = dxdxutrm +dydyutrm +source(i,j)
					rsum = rsum +errchk*errchk
				end if
			end do
		end do
		! Calculate normalized residual.
		res = sqrt(rsum)/rnorm
		
		! If the residual is less than the tolerance, SUCCESS!
		if (res < tol_test .AND. iter>0) then

if(DADI_DEBUG==1) then
	res_ave = res_ave + res
	iter_ave = iter_ave + iter
	write(*,"('dadi: SUCCESS ter=',I4,1X,'res=',G10.5,1X,'iter_ave=',G10.5,1X, &
			'res_ave=',G10.5,1X, 'del_t=',G10.5)") iter, res, iter_ave, res_ave, del_t
end if

			do j = 0, ncy
				do i = 0, ncx
					u_in(i,j) = u(i,j)
				end do
			end do
			goto 100
		end if

		! Determine ratio used to find the time step change.  If tpbot 
		! is zero but tptop is finite, consider this a case of a large 
		! ratio and act accordingly.  DWH does about the same thing 
		! except he does NOT discard the solution if tpbot=0.
		if (tpbot > 0.0) ratio = tptop / tpbot
		if (tpbot < 1e-30) ratio = 1.0

if( DADI_DEBUG==1) then
	open(1,file='.\output\dadi_test.dat',position='append')
		write(1,"('dadi: iter=',I4,1X, 'res=',G10.5, 1X,'tol=',G10.5,1X,'del_t=',G10.5)") iter, res, tol_test, del_t
		write(1,"(5X,'ratio=',G10.5,1X,'tptop=',G10.5,1X,'tpbot=',G10.5)") ratio, tptop, tpbot
	close(1)
end if
	
		! Get next time step.
		if (ratio < 0.02)      then
    		del_t =del_t * 8.000
		else if (ratio < 0.05) then
    		del_t = del_t *4.000
		else if (ratio < 0.10) then
    		del_t = del_t *2.000
		else if (ratio < 0.30) then
    		del_t = del_t *1.000
		else if (ratio < 0.40) then
    		del_t = del_t *0.500
		else if (ratio < 0.60) then
    		del_t = del_t *0.250
		else 
			! Ratio is too large.
			ndiscard = ndiscard + 1
if(DADI_DEBUG==1) then
	write(*,"('iter=',I4,1X,'res=',G10.5,1X,'ndiscard=',I4)") iter, res, ndiscard
end if
			! Check if too many discards
			if (ndiscard > 20) then
!if(DADI_DEBUG==1) then
	write(*, "('too many dischard: ', I4)") ndiscard
!end if
				do j = 0, ncy
					do i = 0, ncx
						u_in(i,j) = u(i,j)
					end do
				end do
				goto 100
			end if
			
			! Discard by replacing u with what we started with.
			do j = 0, ncy
				do i = 0, ncx
					u(i,j) = ustor(i,j)
				end do
			end do
			
			! Reduce del_t.
			del_t = del_t /16.0
		end if
		del_td = 2.*del_t
	end do

	! Fail if used up the maximum iterations.
if (DADI_DEBUG == 1) then
  write(*,"('FAILED iter>=',I4,1X,'res=',G10.5,1X,'tol=',G10.5,1X,'del_t=',G10.5)") &
		itermax, res, tol_test, del_t
end if

	do j = 0, ncy
		do i = 0, ncx
			u_in(i,j) = u(i,j)
		end do
	end do

	write(*,*) iloop, res, sqrt(rsum)

100	continue
	
end subroutine dadisolve




!  Single Peaceman Rachford Douglas pass with Direchlet 0 c boundary
!  conditions for the equation: 
  
!  dtu = dxdxu + dydyu -s, where s is constant in time.  
  
!  The Crank-Nicolson finite difference approximation to the 
!  above equation leads to the fractional step or 
!  ADI equations: 
  
!  u*(i,j)-(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)] 
!  = un(i,j)+(del_t/2dydy)[un(i,j+1)-2un(i,j)+un(i,j-1)] - (del_t/2)s(i,j) 
  
!  un+1(i,j)-(del_t/2dydy)[un+1(i,j+1)-2un+1(i,j)+un+1(i,j-1)] 
!  = u*(i,j)+(del_t/2dxdx)[u*(i+1,j)-2u*(i,j)+u*(i-1,j)] - (del_t/2)s(i,j) 
subroutine adi(uadi, del_t, ncx, ncy)

	use poisson
	implicit none
	
	integer i, j,ip1, im1, jp1, jm1,ncx,ncy, ii, jj
	integer rho_flag, gradient_flag, bdtype_flag
	real a_x(0:ncx), b_x(0:ncx), c_x(0:ncx), r_x(0:ncx), v_x(0:ncx)
	real a_y(0:ncy), b_y(0:ncy), c_y(0:ncy), r_y(0:ncy), v_y(0:ncy)
	real uadi(0:ncx,0:ncy)
	
	real del_t, dthi
	
	dthi = -2.0/del_t
	ustar = uadi
	
	do j = 0, ncy
        b_x = 1.     ! (i,j)的系数
        c_x = 0.     ! (i+1,j)的系数
        a_x = 0.     ! (i-1,j)的系数
        r_x = 0.     ! 右端项
        
        jm1 = j-1
		if (jm1 < 0) jm1 = 0
        jp1 = j+1
		if (jp1 > ncy) jp1 = ncy
		
		do i = 0, ncx
            rho_flag = ro_boundary_flag(i,j)
            gradient_flag = cmprg_gradient(i,j)
            bdtype_flag = cmprg_bdtype(i,j)

			if (rho_flag > 0) then
				a_x(i) = a_xgeom(i,j)
				b_x(i) = -b_xgeom(i,j)
				c_x(i) = c_xgeom(i,j)

				if (abs(a_xgeom(i,j))<=float_error .AND. abs(b_xgeom(i,j)-1.)<=float_error &
					.AND. abs(c_xgeom(i,j))<=float_error) then
					r_x(i) = -dirichletb(i,j)
				else
					b_x(i) = dthi +b_x(i)
					r_x(i) = dthi*uadi(i,j) -source(i,j) &
							-(a_ygeom(i,j)*uadi(i,jm1)-b_ygeom(i,j)*uadi(i,j)+c_ygeom(i,j)*uadi(i,jp1))
				end if
			else
				a_x(i) = 0.
                b_x(i) = 1.
                c_x(i) = 0.
                r_x(i) = uadi(i,j)
			end if
        end do
        
		!Solve tridiagonal system.
		if (x_solver_flag(j,1)==0) then
			call ordinary_tridiag(ncx+1, a_x, b_x, c_x, r_x, v_x)
		else if (x_solver_flag(j,1)==1) then
			call cyclic_tridag(ncx, a_x(0:ncx-1), b_x(0:ncx-1), c_x(0:ncx-1), r_x(0:ncx-1), v_x(0:ncx-1))
			v_x(ncx) = v_x(0)
		else
			write(*,*) "error in X_SOLVER_FLAG"
			stop
		end if

		! Copy solution into ustar. 
		do i = 0, ncx
			ustar(i,j) =v_x(i)
		end do
    end do
	
    do i = 0, ncx-1
        b_y = 1.     ! (i,j)的系数
        c_y = 0.     ! (i,j+1)的系数
        a_y = 0.     ! (i,j-1)的系数
        r_y = 0.     ! 右端项
        
        im1 = i-1
		if (im1 == -1) im1 = ncx-1
        ip1 = i+1
		if (ip1 == ncx) ip1 = 0
        
        do j = 0, ncy
            rho_flag = ro_boundary_flag(i,j)
            gradient_flag = cmprg_gradient(i,j)
            bdtype_flag = cmprg_bdtype(i,j)

			if (rho_flag > 0) then
				a_y(j) = a_ygeom(i,j)
				b_y(j) = -b_ygeom(i,j)
				c_y(j) = c_ygeom(i,j)

				if (abs(a_ygeom(i,j))<=float_error .AND. abs(b_ygeom(i,j)-1.)<=float_error &
					.AND. abs(c_ygeom(i,j))<=float_error) then
					r_y(j) = -dirichletb(i,j)
				else
					b_y(j) = dthi +b_y(j)
					r_y(j) = dthi*ustar(i,j) -source(i,j) &
							-(a_xgeom(i,j)*ustar(im1,j)-b_xgeom(i,j)*ustar(i,j)+c_xgeom(i,j)*ustar(ip1,j))
				end if
			else
				a_y(j) = 0.
                b_y(j) = 1.
                c_y(j) = 0.
                r_y(j) = ustar(i,j)
			end if
        end do
        
		!Solve tridiagonal system.
!		if (y_solver_flag(i,1)==0) then
			call ordinary_tridiag(ncy+1, a_y, b_y, c_y, r_y, v_y)
!		else if (y_solver_flag(i,1)==1) then
!			call cyclic_tridag(ncy, a_y, b_y, c_y, r_y, v_y)
!			v_y(ncy) = v_y(0)
!		else
!			write(*,*) "error in Y_SOLVER_FLAG"
!			stop
!		end if

		! Copy solution into uadi. 
		do j = 0, ncy
			uadi(i,j) = v_y(j)
		end do
	end do

	uadi(ncx,:) = uadi(0,:)

end subroutine adi


!  Tridiagonal field solver:  
  
!  | b0 c0 0                      | | u0 |     | r0 |
!  | a1 b1 c1                     | | u1 |     | r1 |
!  |      ............            | | .  |  =  | .  |
!  |               an-2 bn-2 cn-2 | |un-2|     |rn-2|
!  |               0    an-1 bn-1 | |un-1|     |rn-1|            
subroutine ordinary_tridiag(n,a,b,c,d,x)
	implicit none
    integer :: n
    real :: a(0:n-1),b(0:n-1),c(0:n-1),d(0:n-1),x(0:n-1)
	real, allocatable :: gam(:)
    integer i
    real bet

	allocate(gam(0:n-1))
    
    ! Decomposition and forward substitution.  
    bet = b(0)
    x(0)= d(0)/bet
    
    do i = 1,n-1
        gam(i) = c(i-1)/bet
        bet = b(i) - a(i)*gam(i)
        x(i) = (d(i) -a(i)*x(i-1))/bet
    end do

    !Back substitution.
    do i = n-2,0,-1
	    x(i) =x(i)- gam(i+1)*x(i+1)
    end do

	deallocate(gam)

end subroutine ordinary_tridiag



!   Cyclic tridiagonal field solver using generalized Thomas algorithm
!   (Can also be used to solve the ordinary tridiagonal equation)
!   参考文献: 王兴波, 钟志华, 求解周期性三对角方程组的广义Thomas算法. 计算力学学报, 21(1), 73-76, 2004
  
!  | b0 c0 0                 a0   | | x0 |     | d0 |
!  | a1 b1 c1                     | | x1 |     | d1 |
!  |      ............            | | .  |  =  | .  |
!  |               an-2 bn-2 cn-2 | |xn-2|     |dn-2|
!  | cn-1          0    an-1 bn-1 | |xn-1|     |dn-1|
subroutine cyclic_tridag(n,a,b,c,d,x)
    implicit none
    integer n, i
    real a(0:n-1),b(0:n-1),c(0:n-1),d(0:n-1),x(0:n-1)
    real, allocatable :: q(:), u(:), t(:), h(:), g(:)
    real p
    
    allocate(q(0:n-1))
    allocate(u(0:n-1))
    allocate(t(0:n-1))
    allocate(h(0:n-1))
    allocate(g(0:n-1))
    
    p = b(0)
    q(0)= -c(0)/p
    t(0)= -a(0)/p
    u(0)= d(0)/p
    
    ! 追过程
    do i = 1,n-1
        p = a(i)*q(i-1)+b(i)
	    u(i) = (d(i)-a(i)*u(i-1))/p
	    q(i) = -c(i)/p
	    t(i) = -a(i)*t(i-1)/p
	end do
    
	h(n-1) = (d(n-1)-a(n-1)*u(n-2))/(a(n-1)*(q(n-2)+t(n-2))+b(n-1))
	g(n-1) = -1*c(n-1)/(a(n-1)*(q(n-2)+t(n-2))+b(n-1))
	h(n-2) = u(n-2)+(q(n-2)+t(n-2))*h(n-1)
	g(n-2) = (q(n-2)+t(n-2))*g(n-1)
    
    do i = n-3,1,-1
        h(i) = u(i)+q(i)*h(i+1)+t(i)*h(n-1)
	    g(i) = q(i)*g(i+1)+t(i)*g(n-1)
	end do
	
	! 赶过程
    x(0) = (d(0)-c(0)*h(1)-a(0)*h(n-1))/(b(0)+c(0)*g(1)+a(0)*g(n-1))
    
    do i = 1,n-1
        x(i) = h(i)+g(i)*x(0)
    end do
    
    deallocate(q)
    deallocate(u)
    deallocate(t)
    deallocate(h)
    deallocate(g)
    
end subroutine cyclic_tridag

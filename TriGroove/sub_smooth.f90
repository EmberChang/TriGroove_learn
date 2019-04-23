!***************************************************************!
!* Smoothing theArray using the 1-2-1 method with the proper   *!
!* boundary conditions to conserve charge.                     *!

subroutine smooth(matrix,ncx,ncy)
	use mod_dadi
	use mod_smooth
	implicit none

	integer ncx, ncy
	integer i, j
	integer bdtype_flag, gradient_flag, rho_flag
	real matrix(0:ncx, 0:ncy), temp(0:ncx, 0:ncy)

	if(init_smth_flag) then
		call init_smth_coeff(ncx,ncy) ! 初始化平滑用系数矩阵
	    init_smth_flag = 0
    end if

	temp = 0
	do j = 0, ncy
		do i = 0, ncx
            bdtype_flag = cmprg_bdtype(i,j)
			gradient_flag = cmprg_gradient(i,j)
			rho_flag = ro_boundary_flag(i,j)

			if (bdtype_flag == 3) then
				temp(i,j) = sx_coeff(i,i-1,j)*matrix(i-1,j) +sx_coeff(i,i,j)*matrix(i,j) +sx_coeff(i,i+1,j)*matrix(i+1,j)

			else if (bdtype_flag == 2 .OR. bdtype_flag == 5) then
				if (gradient_flag == 2) then
					temp(i,j) = sx_coeff(i,i,j)*matrix(i,j) +sx_coeff(i,i+1,j)*matrix(i+1,j)
				else if (gradient_flag == 4) then
					temp(i,j) = sx_coeff(i,i-1,j)*matrix(i-1,j) +sx_coeff(i,i,j)*matrix(i,j)
				else if (gradient_flag == 5 .OR. gradient_flag == 6) then
					temp(i,j) = sx_coeff(i,i-1,j)*matrix(i-1,j) +sx_coeff(i,i,j)*matrix(i,j) +sx_coeff(i,i+1,j)*matrix(i+1,j)
				end if
				
			else if (bdtype_flag == 4) then
				if (gradient_flag == 2) then
					if (i /= 0) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						temp(0,j) = sx_coeff(0,0,j)*matrix(0,j) +sx_coeff(0,1,j)*matrix(1,j) +sx_coeff(0,ncx-1,j)*matrix(ncx-1,j)
					end if
				else if (gradient_flag == 4) then
					if (i /= ncx) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						temp(ncx,j) = temp(0,j)
					end if
				else if (gradient_flag == 5 .OR. gradient_flag == 6) then
					temp(i,j) = sx_coeff(i,i-1,j)*matrix(i-1,j) +sx_coeff(i,i,j)*matrix(i,j) +sx_coeff(i,i+1,j)*matrix(i+1,j)
				end if
			
			else if (bdtype_flag == 22) then
				if (rho_flag == 3) then
					temp(i,j) = sx_coeff(i,i-1,j)*matrix(i-1,j) +sx_coeff(i,i,j)*matrix(i,j) +sx_coeff(i,i+1,j)*matrix(i+1,j)
				else if (rho_flag == 4) then
					if (gradient_flag==7 .OR. gradient_flag==8) then
						temp(i,j) = sx_coeff(i,i,j)*matrix(i,j) +sx_coeff(i,i+1,j)*matrix(i+1,j)
					else if (gradient_flag==9 .OR. gradient_flag==10) then
						temp(i,j) = sx_coeff(i,i-1,j)*matrix(i-1,j) +sx_coeff(i,i,j)*matrix(i,j)
					end if
				else
					write (*, "('Error: Exception on rho_flag value ', I8)") rho_flag
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==24 .OR. bdtype_flag==42) then
				if (gradient_flag == 7) then
					if (i /= 0) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						temp(0,j) = sx_coeff(0,0,j)*matrix(0,j) +sx_coeff(0,1,j)*matrix(1,j) +sx_coeff(0,ncx-1,j)*matrix(ncx-1,j)
					end if
				else if (gradient_flag == 9) then
					if (i /= ncx) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						temp(ncx,j) = temp(0,j)
					end if
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==45 .OR. bdtype_flag==54) then
				if (gradient_flag == 8) then
					if (i/=0 .AND. j/=ncy) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						temp(0,j) = sx_coeff(0,0,j)*matrix(0,j) +sx_coeff(0,1,j)*matrix(1,j) +sx_coeff(0,ncx-1,j)*matrix(ncx-1,j)
					end if
				else if (gradient_flag == 10) then
					if (i/=ncx .AND. j/=ncy) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						temp(ncx,j) = temp(0,j)
					end if
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			else if (bdtype_flag == 0) then
				temp(i,j) = 0

			else
				write (*, "('Error: Exception on boundary type ', I8)") bdtype_flag
				stop

            end if
        end do

    end do

	matrix = 0
    do i = 0, ncx
		do j = 0, ncy
            bdtype_flag = cmprg_bdtype(i,j)
			gradient_flag = cmprg_gradient(i,j)
			rho_flag = ro_boundary_flag(i,j)

			if (bdtype_flag == 3) then
				matrix(i,j) = sy_coeff(j,j-1,i)*temp(i,j-1) +sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)

			else if (bdtype_flag == 2 .OR. bdtype_flag == 5) then
				if (gradient_flag == 5) then
					matrix(i,j) = sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)
				else if (gradient_flag == 6) then
					matrix(i,j) = sy_coeff(j,j-1,i)*temp(i,j-1) +sy_coeff(j,j,i)*temp(i,j)
				else if (gradient_flag == 2 .OR. gradient_flag == 4) then
					matrix(i,j) = sy_coeff(j,j-1,i)*temp(i,j-1) +sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)
				end if
				
			else if (bdtype_flag == 4) then
				if (gradient_flag == 5) then
					if (j /= 0) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						matrix(i,0) = sy_coeff(0,0,i)*temp(i,0) +sy_coeff(0,1,i)*temp(i,1) +sy_coeff(0,ncy-1,i)*temp(i,ncy-1)
					end if
				else if (gradient_flag == 6) then
					if (j /= ncy) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						matrix(i,ncy) = matrix(i,0)
					end if
				else if (gradient_flag == 2 .OR. gradient_flag == 4) then
					matrix(i,j) = sy_coeff(j,j-1,i)*temp(i,j-1) +sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)
				end if
			
			else if (bdtype_flag == 22) then
				if (rho_flag == 3) then
					matrix(i,j) = sy_coeff(j,j-1,i)*temp(i,j-1) +sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)
				else if (rho_flag == 4) then
					if (gradient_flag==7 .OR. gradient_flag==9) then
						matrix(i,j) = sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)
					else if (gradient_flag==8 .OR. gradient_flag==10) then
						matrix(i,j) = sy_coeff(j,j-1,i)*temp(i,j-1) +sy_coeff(j,j,i)*temp(i,j)
					end if
				else
					write (*, "('Error: Exception on rho_flag value ', I8)") rho_flag
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==24 .OR. bdtype_flag==42) then
				if (gradient_flag == 7 .OR. gradient_flag == 9) then
					matrix(i,j) = sy_coeff(j,j,i)*temp(i,j) +sy_coeff(j,j+1,i)*temp(i,j+1)
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==45 .OR. bdtype_flag==54) then
				if (gradient_flag == 8) then
					if (i/=0 .AND. j/=ncy) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						matrix(0,ncy) = sy_coeff(ncy,ncy-1,0)*temp(0,ncy-1) +sy_coeff(ncy,ncy,0)*temp(0,ncy)
					end if
				else if (gradient_flag == 10) then
					if (i/=ncx .AND. j/= ncy) then
						write (*, "('Error: Exception of periodic boundary')")
						stop
					else
						matrix(ncx,ncy) = matrix(0,ncy)
					end if
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			else if (bdtype_flag == 0) then
				matrix(i,j) = 0

			else
				write (*, "('Error: Exception on boundary type ', I8)") bdtype_flag
				stop

            end if
        end do

    end do

end subroutine smooth


subroutine init_smth_coeff(ncx,ncy)
	use mod_smooth
	use mod_dadi
	use grid
	implicit none

	integer i, j
	integer gradient_flag, bdtype_flag, bdtype_flag_1, rho_flag
	integer ncx, ncy
	real alpha(0:ncx, 0:ncy)

	sx_coeff = 0.
	sy_coeff = 0.
		
	do j = 0, ncy
		do i = 0, ncx
            bdtype_flag = cmprg_bdtype(i,j)
			gradient_flag = cmprg_gradient(i,j)
			rho_flag = ro_boundary_flag(i,j)

			if (bdtype_flag == 3) then
				sx_coeff(i,i-1,j) = 0.25
				sx_coeff(i,i,j) = 0.5
				sx_coeff(i,i+1,j) = 0.25

			else if (bdtype_flag == 2 .OR. bdtype_flag == 5) then
				if (gradient_flag == 2) then
					sx_coeff(i,i,j) = 0.5
					sx_coeff(i,i+1,j) = 0.5
				else if (gradient_flag == 4) then
					sx_coeff(i,i-1,j) = 0.5
					sx_coeff(i,i,j) = 0.5
				else if (gradient_flag == 5 .OR. gradient_flag == 6) then
					sx_coeff(i,i-1,j) = 0.25
					sx_coeff(i,i,j) = 0.5
					sx_coeff(i,i+1,j) = 0.25
				end if
				
			else if (bdtype_flag == 4) then
				if (gradient_flag == 2) then
					sx_coeff(0,ncx-1,j) = 0.25
					sx_coeff(0,0,j) = 0.5
					sx_coeff(0,1,j) = 0.25
				else if (gradient_flag == 4) then
					continue
				else if (gradient_flag == 5 .OR. gradient_flag == 6) then
					sx_coeff(i,i-1,j) = 0.25
					sx_coeff(i,i,j) = 0.5
					sx_coeff(i,i+1,j) = 0.25
				end if
			
			else if (bdtype_flag==22) then
				if (rho_flag == 3) then
					sx_coeff(i,i-1,j) = 0.25
					sx_coeff(i,i,j) = 0.5
					sx_coeff(i,i+1,j) = 0.25
				else if (rho_flag == 4) then
					if (gradient_flag==7 .OR. gradient_flag==8) then
						sx_coeff(i,i,j) = 0.5
						sx_coeff(i,i+1,j) = 0.5
					else if (gradient_flag==9 .OR. gradient_flag==10) then
						sx_coeff(i,i-1,j) = 0.5
						sx_coeff(i,i,j) = 0.5
					end if
				else
					write (*, "('Error: Exception on rho_flag value ', I8)") rho_flag
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==24 .OR. bdtype_flag==42) then
				if (gradient_flag == 7) then
					sx_coeff(0,ncx-1,j) = 0.25
					sx_coeff(0,0,j) = 0.5
					sx_coeff(0,1,j) = 0.25
				else if (gradient_flag == 9) then
					continue
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==45 .OR. bdtype_flag==54) then
				if (gradient_flag == 8) then
					sx_coeff(0,ncx-1,ncy) = 0.25
					sx_coeff(0,0,ncy) = 0.5
					sx_coeff(0,1,ncy) = 0.25
				else if (gradient_flag == 10) then
					continue
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			else if (bdtype_flag == 0) then
				continue

			else
				write (*, "('Error: Exception on boundary type ', I8)") bdtype_flag
				stop

            end if
        end do

    end do

	do i = 0, ncx
		do j = 0, ncy
            bdtype_flag = cmprg_bdtype(i,j)
			gradient_flag = cmprg_gradient(i,j)
			rho_flag = ro_boundary_flag(i,j)

			if (bdtype_flag == 3) then
				sy_coeff(j,j-1,i) = 0.25
				sy_coeff(j,j,i) = 0.5
				sy_coeff(j,j+1,i) = 0.25

			else if (bdtype_flag == 2 .OR. bdtype_flag == 5) then
				if (gradient_flag == 5) then
					sy_coeff(j,j,i) = 0.5
					sy_coeff(j,j+1,i) = 0.5
				else if (gradient_flag == 6) then
					sy_coeff(j,j-1,i) = 0.5
					sy_coeff(j,j,i) = 0.5
				else if (gradient_flag == 2 .OR. gradient_flag == 4) then
					sy_coeff(j,j-1,i) = 0.25
					sy_coeff(j,j,i) = 0.5
					sy_coeff(j,j+1,i) = 0.25
				end if
				
			else if (bdtype_flag == 4) then
				if (gradient_flag == 5) then
					sy_coeff(0,ncy-1,i) = 0.25
					sy_coeff(0,0,i) = 0.5
					sy_coeff(0,1,i) = 0.25
				else if (gradient_flag == 6) then
					continue
				else if (gradient_flag == 2 .OR. gradient_flag == 4) then
					sy_coeff(j,j-1,i) = 0.25
					sy_coeff(j,j,i) = 0.5
					sy_coeff(j,j+1,i) = 0.25
				end if
			
			else if (bdtype_flag==22) then
				if (rho_flag == 3) then
					sy_coeff(j,j-1,i) = 0.25
					sy_coeff(j,j,i) = 0.5
					sy_coeff(j,j+1,i) = 0.25
				else if (rho_flag == 4) then
					if (gradient_flag==7 .OR. gradient_flag==9) then
						sy_coeff(j,j,i) = 0.5
						sy_coeff(j,j+1,i) = 0.5
					else if (gradient_flag==8 .OR. gradient_flag==10) then
						sy_coeff(j,j-1,i) = 0.5
						sy_coeff(j,j,i) = 0.5
					end if
				else
					write (*, "('Error: Exception on rho_flag value ', I8)") rho_flag
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==24 .OR. bdtype_flag==42) then
				if (gradient_flag == 7 .OR. gradient_flag == 9) then
					sy_coeff(j,j,i) = 0.5
					sy_coeff(j,j+1,i) = 0.5
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			! 这里有经验的因素，在模型边界发生变化时要注意
			else if (bdtype_flag==45 .OR. bdtype_flag==54) then
				if (gradient_flag == 8) then
					sy_coeff(ncy,ncy-1,0) = 0.5
					sy_coeff(ncy,ncy,0) = 0.5
				else if (gradient_flag == 10) then
					continue
				else
					write (*, "('Error: Exception on 内角点边界 with boundary type ', I8)") gradient_flag
					stop
				end if

			else if (bdtype_flag == 0) then
				continue

			else
				write (*, "('Error: Exception on boundary type ', I8)") bdtype_flag
				stop

            end if
        end do

    end do

end subroutine init_smth_coeff
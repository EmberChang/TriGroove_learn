SUBROUTINE GRADIENT(ncx,ncy)

	use constant
	use global
	use poisson
	use wall
	use mod_dadi
	implicit none
	integer ncx, ncy
	integer i, j, k, ipre, jpre

	ez = 0.
	er = 0.
	DO ipre = 0, ncx
		DO jpre = 0, ncy
			if (cmprg_bdtype(ipre,jpre) == 3) then
				ez(ipre,jpre) = (phi(ipre-1,jpre) - phi(ipre+1,jpre)) /(2.*dx)
				er(ipre,jpre) = (phi(ipre,jpre-1) - phi(ipre,jpre+1)) /(2.*dy)
			else if (cmprg_bdtype(ipre,jpre) == 4) then
				er(ipre,jpre) = (phi(ipre,jpre-1) - phi(ipre,jpre+1)) /(2.*dy)
				if (ipre==0 .OR. ipre==ncx) then
					ez(ipre,jpre) = (phi(ncx-1,jpre) - phi(1,jpre)) /(2.*dx)
				else
					write (*, "('Error: Exception on Grid Index')")
					stop
				end if
			end if
		END DO
	END DO

	do i = 1, num_a
		do j = 1, 4
		    if(ib(i,j) == 2) then
				if(j < 3) then
					do k = 0, nr(i)
						ipre = gcb(i,j)
						jpre = gcb(i,3) +k
						if (j == 1) then
							ez(ipre,jpre) = collectq(i,j,k) /2.
							if (cmprg_bdtype(ipre,jpre) == 2) then
								er(ipre,jpre) = (phi(ipre,jpre-1) - phi(ipre,jpre+1)) /(2.*dy)
							else if (cmprg_bdtype(ipre,jpre) == 22) then
								continue
							else
								write (*, "('Error: Exception on BOUNDARY value ', F8.3)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						else
							ez(ipre,jpre) = -collectq(i,j,k) /2.
							if (cmprg_bdtype(ipre,jpre) == 2) then
								er(ipre,jpre) = (phi(ipre,jpre-1) - phi(ipre,jpre+1)) /(2.*dy)
							else if (cmprg_bdtype(ipre,jpre) == 22) then
								continue
							else
								write (*, "('Error: Exception on BOUNDARY value ', F8.3)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						end if
					end do
				else
					do k = 0, nz(i)
						ipre = gcb(i,1) +k
						jpre = gcb(i,j)

						if (j == 3) then
							er(ipre,jpre) = collectq(i,j,k) /2.
							if (cmprg_bdtype(ipre,jpre) == 2) then
								ez(ipre,jpre) = (phi(ipre-1,jpre) - phi(ipre+1,jpre)) /(2.*dx)
							else if (cmprg_bdtype(ipre,jpre) == 22) then
								continue
							else if (cmprg_bdtype(ipre,jpre) == 24 .OR. cmprg_bdtype(ipre,jpre) == 42) then
								if (ipre==0 .OR. ipre==ncx) then
									ez(ipre,jpre) = (phi(ncx-1,jpre) - phi(1,jpre)) /(2.*dx)
								else
									write (*, "('Error: Exception on Grid Index')")
									stop
								end if
							else
								write (*, "('Error: Exception on BOUNDARY value ', F8.3)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						else
							er(ipre,jpre) = -collectq(i,j,k) /2.
							if (cmprg_bdtype(ipre,jpre) == 2) then
								ez(ipre,jpre) = (phi(ipre-1,jpre) - phi(ipre+1,jpre)) /(2.*dx)
							else if (cmprg_bdtype(ipre,jpre) == 22) then
								continue
							else
								write (*, "('Error: Exception on BOUNDARY value ', F8.3)") cmprg_bdtype(ipre,jpre)
								stop
							end if
						end if
					end do
				end if

			else if(ib(i,j) == 5) then
				if(j == 2) then
					write (*, "('Under Construction !!!')")
					stop
				else if(j == 4) then
					do k = 0, nz(i)
						ipre = gcb(i,1) +k
						jpre = gcb(i,j)
						er(ipre,jpre) = source(ipre,jpre)*dy/2. + phi(ipre,jpre-1)/dy ! dirichlet
!						er(ipre,jpre) = 0.	! neumann
						if (cmprg_bdtype(ipre,jpre) == 5) then
							ez(ipre,jpre) = (phi(ipre-1,jpre) - phi(ipre+1,jpre)) /(2.*dx)
						else if (cmprg_bdtype(ipre,jpre) == 45 .OR. cmprg_bdtype(ipre,jpre) == 54) then
							if (ipre==0 .OR. ipre==ncx) then
								ez(ipre,jpre) = (phi(ncx-1,jpre) - phi(1,jpre)) /(2.*dx)
							else
								write (*, "('Error: Exception on Grid Index')")
								stop
							end if
						else
							write (*, "('Error: Exception on BOUNDARY value ', F8.3)") cmprg_bdtype(ipre,jpre)
							stop
						end if
					end do
				else
					write (*, "('Error: Exception on BOUNDARY value')")
					stop
				end if
			end if

		end do
	end do

END SUBROUTINE GRADIENT

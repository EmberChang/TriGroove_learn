SUBROUTINE dielectric(INUM_A, j, IBOUND, K, insp)
	USE constant
	use global
	use wall
	use neutral
	implicit none

	real fractionReflected,fractionScattered,fractioncollect,SECONDARYONE,SECONDARYTWO
	real Te_PRE, VE_SEC, VE_PRE
	real ranum, ranum1
	INTEGER INUM_A, IBOUND, j, K, insp, new

	if (insp==1) then
		
		if(energy_analyze_diag==1) then
			call linear_weighting(insp, inum_a, ibound, k, j, energy_loss, &
									num_a, nsp, 4, imax, e*t(insp,j), 1)
		end if

		call linear_weighting(insp, inum_a, ibound, k, j, IncNum, &
								num_a, nsp, 4, imax, 1., 1)
		
		Te_PRE = t(insp,j)
		VE_PRE = SQRT(vz(insp,j)**2+vr(insp,j)**2+vtheta(j)**2)
		VE_SEC = SQRT(TE_SEC/t_parameter(insp)) ! VELOCITY OF SECONDARY ELECTRON
		
		fractioncollect		 = 0.5*EXP(-Te_PRE**2/43.46**2)
		fractionReflected	 = 0.5*EXP(-Te_PRE**2/30**2)
		secondarytwo		 = 1 -EXP(-Te_PRE**2/127.8958**2)
		secondaryone		 = 1 -fractioncollect -fractionReflected -secondarytwo
		if(secondaryone < 0) secondaryone = 0.

		CALL RANDOM(RANUM)
		If (RANUM < fractioncollect) THEN
		! electron is absorbed by the wall
			call linear_weighting(insp, inum_a, ibound, k, j, collect, &
									num_a, nsp, 4, imax, 1., 1)
			z(insp,j) = -100

		ELSE IF(RANUM >= fractioncollect .AND. RANUM < (fractioncollect+fractionReflected)) THEN
		! 电子被反射
			ICOLLIDE_WALL(insp,j) = 1
			
			call random(ranum1)
			ranum1 = 1.
			ve_pre = ve_pre * ranum1
			IF(IBOUND == 1) THEN
				z(insp,j) = CB(INUM_A, K) + 0.0001*DZ(INUM_A)*(3.-2.*K)
				CALL RVELC(vr(insp,j), vtheta(j), vz(insp,j), VE_PRE, K, 1)
			ELSE
				r(insp,j) = CB(INUM_A, K+2) + 0.0001*DR(INUM_A)*(3.-2.*K)
				CALL RVELC(vz(insp,j), vtheta(j), vr(insp,j), VE_PRE, K, 1)
			END IF
			T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2+Vtheta(j)**2)
			SEE_flag(j) = 1

			if(energy_analyze_diag == 1) then
				call linear_weighting(insp, inum_a, ibound, k, j, energy_loss, &
										num_a, nsp, 4, imax, e*t(insp,j), 2)
			end if

			call linear_weighting(insp, inum_a, ibound, k, j, EmitNum, &
									num_a, nsp, 4, imax, 1., 1)
			
		ELSE IF(RANUM >= (fractioncollect+fractionReflected)) THEN
		! EMIT ONE OR TWO SECONDARY ELECTRON
			ICOLLIDE_WALL(insp,j) = 1
			
			IF(IBOUND == 1) THEN
				z(insp,j) = CB(INUM_A, K) + 0.0001*DZ(INUM_A)*(3.-2.*K)
				CALL RVELC(vr(insp,j), vtheta(j), vz(insp,j), VE_SEC, K, 0)
			ELSE
				r(insp,j) = CB(INUM_A, K+2) + 0.0001*DR(INUM_A)*(3.-2.*K)
				CALL RVELC(vz(insp,j), vtheta(j), vr(insp,j), VE_SEC, K, 0)
			END IF
			T(insp,j) = t_parameter(insp)*(VZ(insp,j)**2+VR(insp,j)**2+Vtheta(j)**2)
			SEE_flag(j) = 1
			
			if(energy_analyze_diag == 1) then
				call linear_weighting(insp, inum_a, ibound, k, j, energy_loss, &
										num_a, nsp, 4, imax, e*t(insp,j), 2)
			end if

			call linear_weighting(insp, inum_a, ibound, k, j, EmitNum, &
									num_a, nsp, 4, imax, 1., 1)
			
			IF(RANUM >= (1-SECONDARYTWO)) THEN
			! TWO SECONDARY EDECTRON, we need add a new electron
			! and on surface there will be a new ion
				particle_add_num(insp) = particle_add_num(insp) + 1
				new = num(insp) + particle_add_num(insp)
				ICOLLIDE_WALL(insp,new) = 1
				AT(insp,new) = AT(insp,j)
				
				IF(IBOUND == 1) THEN
					z(insp,new) = CB(INUM_A, K) + 0.0001*DZ(INUM_A)*(3.-2.*K)
					r(insp,new) = r(insp,j)
					r_ini(insp,new) = r(insp,new)
					CALL RVELC(vr(insp,new), vtheta(new), vz(insp,new), VE_SEC, K, 0)
				ELSE
					r(insp,new) = CB(INUM_A, K+2) + 0.0001*DR(INUM_A)*(3.-2.*K)
					r_ini(insp,new) = r(insp,new)
					z(insp,new) = z(insp,j)
					CALL RVELC(vz(insp,new), vtheta(new), vr(insp,new), VE_SEC, K, 0)
				END IF

				T(insp,new) = t_parameter(insp)*(VZ(insp,new)**2 &
											+ VR(insp,new)**2 +Vtheta(new)**2)
				SEE_flag(new) = 1

				call linear_weighting(insp, inum_a, ibound, k, new, collect, &
										num_a, nsp, 4, imax, 1., 2)

				if(energy_analyze_diag == 1) then
					call linear_weighting(insp, inum_a, ibound, k, j, energy_loss, &
											num_a, nsp, 4, imax, e*t(insp,new), 2)
				end if

				call linear_weighting(insp, inum_a, ibound, k, j, EmitNum, &
										num_a, nsp, 4, imax, 1., 1)

			END IF
		END IF

	else if(insp==2) then

		if(energy_analyze_diag==1) then
			call linear_weighting(insp, inum_a, ibound, k, j, energy_loss, &
									num_a, nsp, 4, imax, e*t(insp,j), 1)
		end if

		call linear_weighting(insp, inum_a, ibound, k, j, collect, &
							num_a, nsp, 4, imax, 1., 1)

		z(insp,j)=-100
		
	end if

END SUBROUTINE dielectric
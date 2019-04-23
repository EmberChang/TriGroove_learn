subroutine mainDSMC
	use constant 
	use global
	use DSMC
	implicit none
	integer :: inspa=3
	if(idsmc==0) then
		call init_DSMC
		idsmc=1
	end if
	call indexm(inspa)
	call sample
	call collmr  
end subroutine mainDSMC

subroutine init_DSMC
	use constant
	use global
	use DSMC
	implicit none
	integer insp,mnn,i,j,inum_a,mc
	real ranum
	real, external :: gam
	insp=3
    SP(1)=2.18E-25 !SP(1) is the molecular mass
    SP(2)=4.939E-10 !SP(2) is the molecular diameter 
    SPM(2)=300
    SPM(3)=0.5
    SPM(4)=1
    SPM(1)=PI*SP(2)**2
    SPM(5)=GAM(2.5-SPM(3))
!SPM(1) 碰撞截面
!SPM(2) the reference temperature 
!SPM(3) the viscosity temperature power law 
!SPM(4) the reciprocal of the VSS scattering parameter
!SPM(5) the gamma function of (5/2-viscosity-temperature power law)
	do i=1,imax1
		do j=1,imax2
			mc=(i-1)*imax2+j
			cc(mc)=pi*((r0+j*drr)**2-(r0+(j-1)*drr)**2)*dzr*lmdd**3
			call random(ranum)
			CCG(2,mc)=ranum
			CCG(1,mc)=SPM(1)*300.*SQRT(11600*t_int(insp)/300.)			
!the maximum value of the (rel. speed)*(cross-section) is set to a reasonable, 
!but low, initial value and will be increased as necessary
		end do
	end do
end subroutine init_DSMC

SUBROUTINE INDEXM(insp)
	use constant
	use global
	use DSMC
	use grid
	implicit none

	integer i,j,k,insp,MC,ipret,jpret
	ic = 0
	j = 1
	do while(j <= num(insp))
		ipret = INT((z(insp,j) - z1) / dzr) +1
		jpret = INT((r(insp,j) - r1) / drr) +1	
		if(ipret>imax1 .OR. ipret<0 .OR. jpret>imax2 .OR. jpret<0) then
		    call remove(insp,j)
		    goto 100
		end if
    	MC = (ipret-1)*imax2 + jpret
	    isc(j) = MC	! 粒子j所在的网格编号是MC
		Ic(2,Mc) = IC(2,Mc) + 1	! IC(2,MC)：编号为MC的网格中的粒子数
	100	j = j + 1
	end do

	j = 0
	do i = 1, mnc
		IC(1,i) = j ! IC(1,i)：第i个网格之前的粒子总数
		j = j + IC(2,i)
	end do

	IC(2,:) = 0

	do j = 1, num(insp)
		Mc = ISC(j)
		IC(2,MC) = IC(2,MC) + 1
		k = IC(1,MC) + IC(2,MC)
		IR(k) = j ! IR：编号为j的粒子根据网格顺序的编号为K
	end do
end subroutine

SUBROUTINE SAMPLE
!--sample the molecules in the flow.\
	use constant
	use global
	use DSMC
	use grid
	implicit none
	integer i,j,k,ipercel,ip,insp
	insp=3
	ncol=0
	sept=0
	do i=1,mnc
		cs(1,i)=1.e-6
		do j=2,7
			cs(j,i)=0.
		end do
	end do
	do i=1,mnc
		ipercel=IC(2,i)
		if(ipercel >0 ) then
			do j=1,ipercel
				k=IC(1,i)+ipercel
				ip=IR(k)
				cs(1,i)=cs(1,i)+1 !--CS()是各个网格中的取样信息
				cs(2,i)=cs(2,i)+vz(insp,ip)
				cs(3,i)=cs(3,i)+vr(insp,ip)
				cs(4,i)=cs(4,i)
				cs(5,i)=cs(5,i)+vz(insp,ip)**2
				cs(6,i)=cs(6,i)+vr(insp,ip)**2
				cs(7,i)=cs(7,i)
			end do
		end if
	end do
end subroutine sample

SUBROUTINE COLLMR !--calculates collisions
	use constant
	use global
	use neutral
	use DSMC
	use grid
	implicit none
	integer i,j,k,nsel,ip1,ip2,isel,insp
	real sn,asel,ranum,vccm(3),vrcp(3),vrc(3),vrr,a,b,c,oc,sc,d,cvr,cvm,avn
	insp=3
	do i=1,mnc
		sn=cs(1,i)
		avn=ic(2,i) !--AVN is the average number of group MM molecules in the cell
		asel=0.5*Ic(2,i)*avn*fnum(1)*fnuma(1)*ccg(1,i)*dt*atomt*neutraltime/omigp/cc(i)+ccg(2,i)
!--ASEL is the number of pairs to be selected, see eqn (11.5)
		nsel=asel		
		ccg(2,i)=asel-nsel
		if(nsel > 0) then
			if(Ic(2,i) < 2) then
				ccg(2,i)=ccg(2,i)+nsel !--if there are insufficient molecules to calculate collisions, the number NSEL is added to the remainer CCG(2,N)
			else
				cvm=ccg(1,i)
				selt=selt+nsel
				do isel=1,nsel
!--需要考察的网格中的碰撞数
					call random(ranum)
1					k=int(ranum*(Ic(2,i)-0.001))+Ic(1,i)+1
					ip1=IR(k)
					call random(ranum)
					k=int(ranum*(Ic(2,i)-0.001))+Ic(1,i)+1
					ip2=IR(k)
					IF(ip1 == ip2) go to 1
					vrc(1)=vz(insp,ip1)-vz(insp,ip2)
					vrc(2)=vr(insp,ip1)-vr(insp,ip2)
					vrc(3)=0 !--VRC(1 to 3) 三个方向的相对速度
					vrr=sqrt(1.5*(vrc(1)**2+vrc(2)**2+vrc(3)**2)) !vrr 相对速度，二维的时候乘以1.5，因为周向为0
					cvr=vrr*veth*SPM(1)*((2.*kb*SPM(2)/(0.5*SP(1)*VRR))**(SPM(3)-0.5))/SPM(5) !--the collision cross-section is based on eqn (4.63)
					if (cvr > cvm) cvm=cvr !--if necessary, the maximum product in CVM is upgraded
					call random(ranum)
					if(ranum < cvr/ccg(1,i)) then !--the collision is accepted with the probability of eqn (11.6)
						ncol=ncol+1 !--NCOL是实际发生的碰撞数
						sept=sept+sqrt((z(insp,ip1)-z(insp,ip2))**2+(r(insp,ip1)-r(insp,ip2))**2)
						vccm(1)=0.5*(vz(insp,ip1)+vz(insp,ip2)) !--VCCM defines the components of the centre-of-mass velocity, eqn (2.1)
						vccm(2)=0.5*(vr(insp,ip1)+vr(insp,ip2))
						vccm(3)=0 !不考虑第周向原子速度
						if(abs(spm(4)-1.) < 1.e-3) then !--use the VHS logic
							call random(ranum)
							B=2*ranum-1 !B is the cosine of a random elevation angle
							a=sqrt(1-b*b)
							vrcp(1)=b*vrr
							call random(ranum)
							c=2*pi*ranum !C is a random azimuth angle
							vrcp(2)=a*cos(c)*vrr
							vrcp(3)=a*sin(c)*vrr
						else !use the VSS logic
							call random(ranum)
							b=2*(ranum**spm(4))-1 !--B is the cosine of the deflection angle for the VSS model, eqn (11.8)
							a=sqrt(1-b*b)
							call random(ranum)
							c=twopi*ranum
							oc=cos(c)
							sc=sin(c)
							d=sqrt(vrc(2)**2+vrc(3)**2)
							if(d>1.e-6) then ! VRCP(1 to 3) are the components of the post-collision relative vel.
						          VRCP(1)=B*VRC(1)+A*SC*D
								  VRCP(2)=B*VRC(2)+A*(vrr*VRC(3)*OC-VRC(1)*VRC(2)*SC)/D
								  VRCP(3)=B*VRC(3)-A*(vrr*VRC(2)*OC+VRC(1)*VRC(3)*SC)/D
							else
						          VRCP(1)=B*VRC(1)
						          VRCP(2)=A*OC*VRC(1)
						          VRCP(3)=A*SC*VRC(1)
							end if
!--the post-collision rel. velocity components are based on eqn (2.22)
						end if
						vz(insp,ip1)=vccm(1)+0.5*vrcp(1)
						vr(insp,ip1)=vccm(2)+0.5*vrcp(2)
						vz(insp,ip2)=vccm(1)-0.5*vrcp(1)
						vr(insp,ip2)=vccm(2)-0.5*vrcp(2)
					end if
				end do
				ccg(1,i)=cvm
			end if
		end if
	end do
end subroutine COLLMR

real FUNCTION GAM(X)
	real x
!--calculates the Gamma function of X.
      A=1.
      Y=X
      IF (Y.LT.1.) THEN
        A=A/Y
      ELSE
50      Y=Y-1
        IF (Y.GE.1.) THEN
          A=A*Y
          GO TO 50
        END IF
      END IF
      GAM=A*(1.-0.5748646*Y+0.9512363*Y**2-0.6998588*Y**3+ &
         0.4245549*Y**4-0.1010678*Y**5)
END FUNCTION gam
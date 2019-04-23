MODULE constant
    use ifport
	integer, parameter :: NUM_A = 22 ! nuber of the region
	integer, parameter :: IMAX = 224, imax1 = 224, imax2 = 210 ! imax为z和r边网格的最大值，imax1为z边网格最大值，imax2为r边网格最大值

!	integer, parameter :: NUM_A = 1 ! nuber of the region ! plane
	integer, parameter :: nsp = 2 ! 带电粒子的种类
	integer, parameter :: NumPar_MAX = 5E6	! MAXIMUM number of each particle species
	integer, parameter :: free_path_diag = 0, trajectory_diag = 0, conductivity_diag = 0, performance_diag = 0, &
							see_diag = 0, energy_analyze_diag = 0, neutral_diag = 0, flux_diag = 0, &
							ionization_rate_diag = 0, hall_diag = 0
	integer, parameter :: case_output = 2.e3
	integer, parameter :: animation = 1 ! 是否存储历史数据文件，用于制作动态过程图像
	integer, parameter :: iatom = 1 ! 是否考虑电离对原子的影响和离子碰撞后变成原子
	integer, parameter :: igas = 1 ! 1 for xenon, 2 for krypton
	integer, parameter :: atomt = 100 ! 离子、原子时间步长与电子时间步长的比值
	real, parameter :: PI = 3.14159265, TWOPI = 6.2831853
	real, parameter :: kb = 1.38E-23 ! 波尔兹曼常数
	real, parameter :: EPSILON0 = 8.85E-12 ! vaccum permitivity
	real, parameter :: e = 1.6022E-19
	real, parameter :: sfactor = 1.
	real, parameter :: mfactor(nsp) = (/1., 1.E-2/)
	
	real, parameter :: float_error = 0.0001 ! 浮点数逻辑判断的误差限
	real, parameter :: co_stability1 = 1., co_stability2 = 0.2
	real, parameter :: NERO = 1E18 ! number density of initial plasma density

	real :: pfactor = 1. ! 真空介电常数不断调节，以保证稳定性
	real :: dt = 1. ! time step, normalized value
END MODULE constant

module global
	use constant

	integer, parameter :: w_groove = 14								! 沟槽半宽，德拜长度标幺化
	integer, parameter :: h_groove = 21								! 沟槽深度，德拜长度标幺化
	integer, parameter :: h_plasma = 84								! 入射边界到平台的高度，德拜长度标幺化
	integer, parameter :: w_wall = 42								! 一侧平台宽度，德拜长度标幺化
	integer, parameter :: lx = 2									! 阶梯的宽度
	integer, parameter :: ly = 3									! 阶梯的高度

	integer, parameter :: cells_per_lmdd = 2
	
	real, parameter :: Q(nsp) = (/-1.6022E-19, 1.6022E-19/)
	real, parameter :: m(nsp) = (/9.1094E-31*mfactor(1), 2.19E-25*mfactor(2)/) ! 电子、离子质量
	real, parameter :: Q1(nsp) = (/-1., 1./)
	real, parameter :: M1(nsp) = (/1., m(2)/m(1)/)
	real, parameter :: half_q_m(nsp) = (/0.5*q1(1)/m1(1), 0.5*q1(2)/m1(2)/)
	real, parameter :: T_int(nsp) = (/10.,1./)
	real, parameter :: T_enter(nsp) = (/10.,1./)
	real, parameter :: te_sec = T_int(1)*0.3 ! ENERGY OF SECONDARY ELECTRON
!	real, parameter :: lmdd0 = 7.43E3*(T_int(1)/NERO)**0.5 ! 初始电子密度变化之后，这条语句也要作相应的修改
	real, parameter :: lmdd0 = (epsilon0*kb*11600.*t_int(1)/(e*e*nero))**0.5 ! Debye length

	REAL :: CB(NUM_A,4) = 0.
!	real :: CB(NUM_A,4) = (/0., pingtai*lmdd0, (pingtai+lx)*lmdd0, pingtai*lmdd0, &
!								pingtai*lmdd0, (pingtai+lx)*lmdd0, (2*pingtai+lx)*lmdd0, (pingtai+lx)*lmdd0, &
!								ly*lmdd0, ly*lmdd0, ly*lmdd0, 0., &
!								(qiaoceng+ly)*lmdd0, (qiaoceng+ly)*lmdd0, (qiaoceng+ly)*lmdd0, ly*lmdd0/)
!	real :: CB(NUM_A,4) = (/0., pingtai*lmdd0, 0., qiaoceng*lmdd0/)  ! plane

	integer :: IB(NUM_A,4) = 0
!	integer, parameter :: IB(NUM_A,4) = (/4, 3, 3, 2, &
!										3, 3, 4, 2, &
!										2, 3, 2, 2, &
!										5, 5, 5, 3/) ! 边界条件
!	integer, parameter :: IB(NUM_A,4) = (/4, 4, 2, 5/) ! 边界条件  ! plane
	! 1 vaccum | 2 dielectric surface | 3 interior | 4 periodic | 5 inlet
	! 6 mirror reflection | 7 anode | 8 dielectric but fix phi
	! 9 轴对称边界条件（粒子弹性反射，电势梯度为0）| 10 无穷远边界条件（粒子消失，电势梯度为0）
	! 11导体边界（粒子消失，电势给定或通过电容公式计算得到）
	
	integer :: NZ(NUM_A) = 0, NR(NUM_A) = 0
!	integer, parameter :: NZ(NUM_A) = (/cells_per_lmdd*pingtai, cells_per_lmdd*lx, &
!										cells_per_lmdd*pingtai, cells_per_lmdd*lx/)
!	integer, parameter :: NZ(NUM_A) = (/cells_per_lmdd*pingtai/) ! plane
!	integer, parameter :: NR(NUM_A) = (/cells_per_lmdd*qiaoceng, cells_per_lmdd*qiaoceng, &
!										cells_per_lmdd*qiaoceng, cells_per_lmdd*ly/)
!	integer, parameter :: NR(NUM_A) = (/cells_per_lmdd*qiaoceng/) ! plane
	
	integer :: gcb(num_a,4) = 0
!	integer, parameter :: gcb(num_a,4) = (/0, cells_per_lmdd*pingtai, cells_per_lmdd*(pingtai+lx), cells_per_lmdd*pingtai, &
!			cells_per_lmdd*pingtai, cells_per_lmdd*(pingtai+lx), cells_per_lmdd*(2*pingtai+lx), cells_per_lmdd*(pingtai+lx), &
!			cells_per_lmdd*ly, cells_per_lmdd*ly, cells_per_lmdd*ly, 0, &
!			cells_per_lmdd*(ly+qiaoceng), cells_per_lmdd*(ly+qiaoceng), cells_per_lmdd*(ly+qiaoceng), cells_per_lmdd*ly/)
!	integer, parameter :: gcb(num_a,4) = (/0, cells_per_lmdd*pingtai, 0, cells_per_lmdd*qiaoceng/) ! plane
	
	integer :: nzt, nrt
!	integer, parameter :: nzt = (2*pingtai+lx)*cells_per_lmdd, nrt = (ly+qiaoceng)*cells_per_lmdd
!	integer, parameter :: nzt = cells_per_lmdd*pingtai, nrt = cells_per_lmdd*qiaoceng ! plane
	
	real :: z1, z2, r0, r1, r2
!	real :: z1 = 0., z2 = (2*pingtai+lx)*lmdd0
!	real :: z1 = 0., z2 = pingtai*lmdd0 ! plane
!	real :: r0 = 0.021-ly*lmdd0, r1=0., r2 = (ly+qiaoceng)*lmdd0
!	real :: r0 = 0.021, r1=0., r2 = qiaoceng*lmdd0 ! plane
	
	real dz(num_a), dr(num_a) ! size of the cell
	real dzr, drr ! In simulation, the cell is uniform, so dzr=dz(i), drr=dr(i)

	real FNUM(num_a) ! 每个区域的粒子权重
	real OMIGP, LMDD, VETH ! PLASMA FREQUENCY，DEBYE LENGTH,THERMOLIZED VELOCITY OF ELECTRON

	! 粒子参数
	real T(nsp, NumPar_MAX) ! 每个粒子的动能
	real vr(nsp, NumPar_MAX), vz(nsp, NumPar_MAX), vtheta(NumPar_MAX)
	real r(nsp, NumPar_MAX), z(nsp, NumPar_MAX)
	real r_ini(nsp,numpar_max)
	real AT(nsp, NumPar_MAX)
	integer ICOLLIDE_WALL(nsp, NumPar_MAX)
	integer SEE_flag(NumPar_MAX)

	! 场参数
	real :: ez(0:IMAX1,0:IMAX2) = 0., er(0:IMAX1,0:IMAX2) = 0.
	
	real :: t_parameter(nsp) = 0. ! 用于计算电子温度时所用的参数，减少计算量

	integer :: num(nsp) = 0, ILOOP = 0, iloopa = 4.e4 ! 每种带电粒子的个数, 循环次数, 计算原子初始分布所需要的循环次数
	integer state ! 1 for neutral computation, 2 for plasma computation
	integer :: ifield1 = 0, ifield2 = 0 ! 用于output子程序，表示输出的是第几个field
end module global

module mcc_para
	use constant
	integer, parameter :: energy_max = 500
	integer :: mccflag = 0
	integer :: ionized_particle = 0 ! 每次电离产生电子和离子数（两者相等）
	real, parameter :: co_bohm = 1./200., co_bohm1 = 1./64.
	real sels(energy_max), sexc(energy_max), sion(energy_max), SV_MAX_E
	real sigmaElastic_kr(energy_max), sigmaExc_kr(energy_max), sigmaIz_kr(energy_max)
	real sece(0:IMAX1,0:IMAX2) ! 电离率，电离产生的电子计数
end module mcc_para

module ion
	use constant
	real vir_ion0
	real fluxqi(NUM_A,0:IMAX) ! ion charge flux
end module ion

module grid
	use constant
	real volume_cell(0:IMAX) ! volume of the cell
	integer, parameter :: mnc = imax1*imax2 ! 计算区域内的网格数
	integer :: IC(2,imax1*imax2) = 0, IR(NumPar_MAX) = 0, ISC(NumPar_MAX) = 0
end module grid
 
module grid_para
	use constant

	real ro_boundary(0:IMAX1,0:IMAX2) ! 用于计算边界密度时的系数
	real :: rho(nsp,0:imax1,0:imax2) = 0.

	!!!!!!!just for output!!!!!!
	real :: rho_g(nsp,0:imax1,0:imax2) = 0.
	real :: vz_g(nsp,0:IMAX1,0:IMAX2) = 0., vr_g(nsp,0:IMAX1,0:IMAX2) = 0., vtheta_g(nsp,0:imax1,0:imax2) = 0.
	real :: egy_g(nsp,0:imax1,0:imax2) = 0., t_g(nsp,0:imax1,0:imax2) = 0.
	real :: vt_g(nsp,0:IMAX1,0:IMAX2) = 0., vt2_g(nsp,0:IMAX1,0:IMAX2) = 0.
	real :: ez_g(0:IMAX1,0:IMAX2) = 0., er_g(0:IMAX1,0:IMAX2) = 0.
	real :: phi_g(0:IMAX1,0:IMAX2) = 0.
	
	real :: rho_g_SEE(nsp,0:imax1,0:imax2) = 0.
	real :: vtheta_g_SEE(nsp,0:imax1,0:imax2) = 0.

	real :: egy_g_normal(nsp,0:imax1,0:imax2) = 0.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!积累了case_output计算步长之后的累计值!!!!!!
	real :: rho_g1(nsp,0:imax1,0:imax2) = 0.
	real :: rhovz_g1(nsp,0:imax1,0:imax2) = 0., rhovr_g1(nsp,0:imax1,0:imax2) = 0., rhovtheta_g1(nsp,0:imax1,0:imax2) = 0.
	real :: rhoegy_g1(nsp,0:imax1,0:imax2) = 0., rhovt_g1(nsp,0:IMAX1,0:IMAX2) = 0., rhovt2_g1(nsp,0:IMAX1,0:IMAX2) = 0.
	real :: ez_g1(0:IMAX1,0:IMAX2)=0., er_g1(0:IMAX1,0:IMAX2)=0.
	real :: phi_g1(0:IMAX1,0:IMAX2) = 0.
	
	real :: rho_g1_SEE(nsp,0:imax1,0:imax2) = 0.
	real :: rhovtheta_g1_SEE(nsp,0:imax1,0:imax2) = 0.

	real :: rhoegy_g1_normal(nsp,0:IMAX1,0:IMAX2) = 0.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module grid_para

module magnetic
	use constant
	real bzt(0:IMAX1,0:IMAX2), brt(0:IMAX1,0:IMAX2), bmax, OMIGC(nsp), octp
end module magnetic

module wall
	use constant
	real, parameter :: t_wall = 700.
	real COLLECT(nsp, NUM_A, 4, 0:IMAX), COLLECTQ(NUM_A, 4, 0:IMAX) ! surface number of each species and net charge density
	real :: IncNum(nsp, NUM_A, 4, 0:IMAX) = 0. ! 入射到壁面的电子数
	real :: EmitNum(nsp, NUM_A, 4, 0:IMAX) = 0. ! 二次电子数
	real :: num_inr(nsp) = 0.
	real :: energy_loss(nsp, num_a, 4, 0:imax)=0.
	real :: conductor_charge = 0. ! 导体边界积累的电荷
	integer :: particle_add_num(nsp) = 0 ! 由于粒子和壁面的相互作用而需要删除和增加的粒子
	integer :: nout(nsp, num_a, 4) = 0 ! number of particle go out of the region per loop
end module wall

module neutral
	use constant
	real :: fnuma(NumPar_MAX) = 500. ! a atom will represent 500 ions and electrons
	real :: fnuma_initial = 2000. ! 原子权重与等离子权重相比的倍数，这个值保持不变，fnuma可能会变
	real, parameter :: atom_flux = 3.e-6*sfactor
	integer, parameter :: neutraltime = 2.e1
end module neutral

module diag
	use constant
	integer, parameter :: ipc_num = 99, free_path_num = 1.e4
	integer :: ipc(ipc_num) = 0 ! a parameter for whether particle 1 collide with boundary, 1 for collide with
	integer :: ir1(ipc_num) = 0, particle_diag = 0 !记录跟踪的粒子在粒子序列中的记号
	integer :: collision_time(free_path_num)=0
	real :: init_z(free_path_num) = 0., init_r(free_path_num) = 0. ! 统计最多一万个电子碰撞的初始位置和经历的时间
	
end module diag

module discharge
	use constant
	real, parameter :: Id=3.*sfactor, Ud=0. ! discharge current and voltage
end module discharge

module DSMC
	use constant
	integer ncol
	integer :: idsmc = 0
	real CS(7, imax1*imax2), spm(5), ccg(2, imax1*imax2)
	real cc(imax1*imax2), selt, sept, sp(2)
end module DSMC

module mod_dadi
	use constant
	real, parameter :: capacitance = 1.e-7
	real PHI(0:IMAX1, 0:IMAX2) ! electric potential
	real phi0(num_a, 4, 0:imax) ! phi on boundary for Direchlet boundary condition
	integer cmprg_gradient(0:imax1, 0:imax2)
	integer cmprg_bdtype(0:imax1, 0:imax2), ro_boundary_flag(0:imax1, 0:imax2)
	integer x_solver_flag(0:imax,1), y_solver_flag(0:imax,1) ! 0, 1 means using ordinary, cyclic tridiagnol solver, respectively.
end module mod_dadi

module mod_smooth
	use constant
	real :: sx_coeff(0:imax1, 0:imax1, 0:imax2) = 0.
	real :: sy_coeff(0:imax2, 0:imax2, 0:imax1) = 0.
	integer :: init_smth_flag = 1
end module mod_smooth

module FermanSEEMod
!	Emitted angular spectrum
	real, parameter :: FSM_alpha = 1.

!	Backscattered electrons
	real, parameter :: FSM_P1e_inf = 0.028
	real, parameter :: FSM_P1e_cusp = 0.45
	real, parameter :: FSM_Ee_cusp = 0.
	real, parameter :: FSM_W = 60.86
	real, parameter :: FSM_p = 1.
	real, parameter :: FSM_sige = 2.
	real, parameter :: FSM_e1 = 0.26
	real, parameter :: FSM_e2 = 2.

!	Rediffused electrons
	real, parameter :: FSM_P1r_inf = 0.276
	real, parameter :: FSM_Er = 0.041
	real, parameter :: FSM_r = 0.104
	real, parameter :: FSM_q = 0.5
	real, parameter :: FSM_r1 = 0.26
	real, parameter :: FSM_r2 = 2.

!	True secondary electrons
	integer, parameter :: FSM_M = 10
	real, parameter :: FSM_dltts_cusp = 2.6
	real, parameter :: FSM_Ets_cusp = 600.
	real, parameter :: FSM_s = 1.54
	real, parameter :: FSM_t1 = 0.66
	real, parameter :: FSM_t2 = 0.8
	real, parameter :: FSM_t3 = 0.7
	real, parameter :: FSM_t4 = 1.
	real, parameter :: FSM_pn(FSM_M) = (/2.5,3.3,2.5,2.5,2.8,1.3,1.5,1.5,1.5,1.5/)
	real, parameter :: FSM_epsin(FSM_M) = (/1.5,1.75,1.,3.75,8.5,11.5,2.5,3.,2.5,3./)

!	Total SEY
	real, parameter :: FSM_dltt_cusp = 2.9

	integer :: iseed = 1	

end module FermanSEEMod

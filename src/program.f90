program sph

	! Main program of project

	!Created by Jannik Zuern on 05/16/2016
	!Last modified: 05/24/2016


	use gnufor2       			!Module: gnuplot fortran visualizations
	! use gnuplot_module  			!better Gnuplot visualization options
	! use gnuplot_module_data		!needed by gnuplot_module
	! use datatypes, only : i4b !needed by gnuplot_module


	use integrate     !Module: definition of integration functions (leapfrog)
	use sphfunctions  !Module: implementation of sph mechanics
	use linkedlists   !Module: Linked lists mechanics and bookkeeping
	use util          !Module: utiliary functions

	implicit none


! ACTUAL PROGRAM STARTING
!-------------------------------------------------------------
!-------------------------------------------------------------

! Define variables
! type(gnuplot_ctrl), pointer 						:: ptr_gctrl           ! gnuplot control pointer
! character(len=1) 												:: debug ! for gnuplotfortran
! integer(i4b) 														:: status=0,ierror=0

integer																	:: nFrames 							! number of frames to calculate
integer																	:: nSteps_per_frame     ! number of steps between shown frame
integer 																:: i,j,k 								! loop iteration variables
integer, allocatable, dimension(:)  	  :: ll 									! linked list array
integer, allocatable, dimension(:,:)  	:: lc 									! linked cell array
double precision, dimension (9)    	   	:: simulation_parameter ! Initialize sim param vector
type(systemstate)                 	    :: sstate								! Simulation state
type(sim_parameter)											:: params               ! parameter of simulation
double precision												:: dt 									! (constant) time step for numerical integration


! Parse input parameter file
call parse_input(simulation_parameter)

! Write contents of parameter array into sim_parameter type (handy usage)
params = initialize_parameters(simulation_parameter)

! initialize particles
call init_particles(sstate,params,ll,lc)

! initialize linked lists
call init_ll(sstate,ll)

! initialize linked cells lists
call init_lc(sstate,params,lc)

! setting up neighbor lists based on placed particles
call setup_neighbour_list(sstate, params, ll,lc)
! call print_neighour_list(sstate, params, ll,lc) ! for debugging

! First integration
print *, "Calculating Step ", 0
call compute_accel(sstate, params, ll,lc)
call leapfrog_start(sstate,params%dt)    ! for first iteration, we must use different leapfrog algorithm
																	       ! since we de not have any previous time step yet

!Invoking gnuplot_init ...
! ptr_gctrl=>gnuplot_init('-persist')
! if(.not.associated(ptr_gctrl)) stop 'Failed to initiate a gnuplot session in program'
! !
! status=gnuplot_setrange(ptr_gctrl,'x',0.0d0,1.0d0)
! status=gnuplot_setrange(ptr_gctrl,'y',0.0d0,1.0d0)
! status=gnuplot_settitle(ptr_gctrl,'Testbench plot')
! status=gnuplot_setaxislabel(ptr_gctrl,'x','x')
! status=gnuplot_setaxislabel(ptr_gctrl,'y','sin,cos,atan,sqrt')
! status=gnuplot_setscale(ptr_gctrl,'x','NLG')
! status=gnuplot_setscale(ptr_gctrl,'y','NLG')
! status=gnuplot_setstyle(ptr_gctrl,'linespoints')
!
! status=gnuplot_plot2d(ptr_gctrl,sstate%nParticles,sstate%x(1:2*sstate%nParticles:2),sstate%x(2:2*sstate%nParticles:2),'sin(x)')


! Simulation loop
nFrames 						= params%nFrames
nSteps_per_frame 		= params%nSteps_per_frame

do i = 1,nFrames
	print *, "Calculating Step ", i, " of " , nFrames

	do j = 1,nSteps_per_frame
		call compute_accel(sstate, params,ll,lc)  !update values for accellerations
		call leapfrog_step(sstate, params%dt) 		!update velocities and positions based on previously calculated accelleration
	end do

	call plot_data_immediately(sstate,i)
	! call write_data_to_file(sstate,i)

end do

! do i = 1,nFrames
! 	call plot_data_from_file(sstate,i)
! end do


! Cleanup
call free_state(sstate)

! don't end program
read *,

end program sph

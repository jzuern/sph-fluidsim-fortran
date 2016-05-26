program sph

	! Main program of project

	!Created by Jannik Zuern on 05/16/2016
	!Last modified: 05/24/2016

	use gnufor2       !Module: gnuplot fortran visualizations
	use integrate     !Module: definition of integration functions (leapfrog)
	use sphfunctions  !Module: implementation of sph mechanics
	use linkedlists   !Module: Linked lists mechanics and bookkeeping
	use util          !Module: utiliary functions

	implicit none





! ACTUAL PROGRAM STARTING
!-------------------------------------------------------------
!-------------------------------------------------------------

! Define variables
integer																	:: nFrames 							! number of frames to calculate
integer																	:: nSteps_per_frame     ! number of steps between shown frame
integer 																:: i,j,k 								! loop iteration variables
integer, allocatable, dimension(:)  	  :: ll 									! linked list array
integer, allocatable, dimension(:,:)  	:: lc 									! linked cell array
double precision, dimension (9)    	   	:: simulation_parameter ! Initialize sim param vector
type(systemstate)                 	    :: sstate								! Simulation state
type(sim_parameter)											:: params               ! parameter of simulation
double precision												:: dt 									! (constant) time step for numerical integration


! wipe data directory from old files



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
call print_neighour_list(sstate, params, ll,lc)

! First integration
print *, "Calculating Step ", 0
call compute_accel(sstate, params, ll,lc)
call leapfrog_start(sstate,params%dt)    ! for first iteration, we must use different leapfrog algorithm
																	! since we de not have any previous time step yet
! call check_state(sstate)   				! currently not working


! Simulation loop
nFrames 						= params%nFrames
nSteps_per_frame 		= params%nSteps_per_frame

do i = 1,nFrames
	print *, "Calculating Step ", i, " of " , nFrames
	do j = 1,nSteps_per_frame
		! print *, "                    Calculating Sub-Step ", j
		call compute_accel(sstate, params,ll,lc) !update values for accellerations
		call leapfrog_step(sstate, params%dt) 												 !update velocities and positions based on previously calculated accelleration
		call check_state(sstate);  												 	   !not working
	end do

	call plot_data_immediately(sstate,i)
	! call write_data_to_file(sstate,i)


end do

do i = 1,nFrames
	! call plot_data_from_file(sstate,i)
end do



! Cleanup
call free_state(sstate)

! don't end program
read *,

end program sph

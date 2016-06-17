program sph

	! Main program of project

	!Created by Jannik Zuern on 05/16/2016
	!Last modified: 06/19/2016


	use gnufor2       			!Module: gnuplot fortran visualizations

	use integrate     !Module: definition of integration functions (leapfrog)
	use sphfunctions  !Module: implementation of sph mechanics
	use linkedlists   !Module: Linked lists mechanics and bookkeeping
	use util          !Module: utiliary functions

	implicit none


! ACTUAL PROGRAM STARTING
!-------------------------------------------------------------
!-------------------------------------------------------------




integer																	:: nFrames 							! number of frames to calculate
integer																	:: nSteps_per_frame     ! number of steps between shown frame
integer 																:: i,j  								! loop iteration variables
integer, allocatable, dimension(:)  	  :: ll 									! linked list array
integer, allocatable, dimension(:,:)  	:: lc 									! linked cell array
type(systemstate)                 	    :: sstate								! Simulation state
type(sim_parameter)											:: params               ! parameter of simulation
integer 																:: t1,t2, clock_rate, clock_max ! for time counting

! start counting of program runtime with intrinsic subroutine system_clock
call system_clock ( t1, clock_rate, clock_max )


! Write contents of parameter array into sim_parameter type (handy usage)
call initialize_parameters(params)

! initialize all particles
call init_particles(sstate,params)


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
call leapfrog_start(sstate,params)    ! for first iteration, we must use different leapfrog algorithm
																	       ! since we de not have any previous time step yet



! Simulation loop
nFrames 						= params%nFrames
nSteps_per_frame 		= params%nSteps_per_frame


! do i = 1,nFrames
! 	print *, "Calculating Step ", i, " of " , nFrames
!
! 	do j = 1,nSteps_per_frame
! 		call compute_accel(sstate, params,ll,lc)  !update values for accellerations
! 		call leapfrog_step(sstate, params) 	    	!update velocities and positions based on previously calculated accelleration
! 	end do
!
!
!   ! we can either plot data immediately (call plot_data_immediately(...) ) ...
! 	! call plot_data_immediately(sstate,i)
!
! 	!.. or we can first write the data to files and plot them later (in separate loop)
! 	call write_data_to_file(sstate,i)
!
! end do


!print out elapsed time:
call system_clock ( t2, clock_rate, clock_max )
write ( *, * ) 'Elapsed time = ', real (t2-t1)/real(clock_rate), 'seconds, which is an average of ' &
														 & ,  real (t2-t1)/real(clock_rate)/real(nFrames), ' seconds per frame.'

!interrupt execution until keyboard interaction
print *, "Press any key to start plotting"
read(*,*)

! this is the plotting loop (if some data exists in folder)
do i = 1,nFrames
	call plot_data_from_file(sstate,i)
end do


! Cleanup
call free_state(sstate)

! end of program
end program sph

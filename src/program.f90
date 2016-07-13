program sph

	! Main program of project

	!Created by Jannik Zuern on 05/16/2016
	!Last modified: 07/05/2016


	use integrate     !Module: definition of integration functions (leapfrog)
	use sphfunctions  !Module: implementation of sph mechanics
	use linkedlists   !Module: Linked lists mechanics and bookkeeping
	use util          !Module: utiliary functions

	use gnuplot_module_data, only : gnuplot_ctrl !gnuplotfortran visualizations

	implicit none


! PROGRAM STARTING
!-------------------------------------------------------------
!-------------------------------------------------------------




integer																	:: nFrames 							! number of frames to calculate
integer																	:: nSteps_per_frame     ! number of steps between shown frame
integer 																:: i,j  								! loop iteration variables
integer, allocatable, dimension(:)  	  :: ll 									! linked list array
integer, allocatable, dimension(:,:,:)  	:: lc 									! linked cell array
type(systemstate)                 	    :: sstate								! Simulation state
type(sim_parameter)											:: params               ! parameter of simulation
integer 																:: t1,t2, clock_rate, clock_max ! for time counting

type(gnuplot_ctrl), pointer 						:: ptr_gctrl ! pointer to gnuplot control unit


! start counting of program runtime with intrinsic subroutine system_clock
call system_clock ( t1, clock_rate, clock_max )

! Write contents of parameter array into sim_parameter type
call initialize_parameters(params)

! initialize all particles
call init_particles(sstate,params)

! initialize linked lists
call init_ll(sstate,ll)

! initialize linked cells lists
call init_lc(sstate,params,lc)

! set up neighbor lists based on placed particles
call setup_neighbour_list(sstate, params, ll,lc)

call write_data_to_file(sstate,0)

! First time-integration
!
! print *, sstate%a


print *, "Calculating frame ", 0
call compute_accel(sstate, params, ll,lc)



call leapfrog_start(sstate,params)    ! for first iteration, we must use different leapfrog algorithm
																	    ! as we de not have any previous time step yet


! call write_data_to_file(sstate,0)

!invoke gnuplot session
! call invoke_gnuplot(ptr_gctrl)


nFrames 						= params%nFrames
nSteps_per_frame 		= params%nSteps_per_frame


! Loop through frames
do i = 1,nFrames
	print *, "Calculating frame ", i, " of " , nFrames

	!Loop through steps for each frame
	do j = 1,nSteps_per_frame
		call compute_accel(sstate, params,ll,lc)  !update values for accellerations
		call leapfrog_step(sstate, params) 	    	!update velocities and positions based on previously calculated accelleration
	end do


  ! Plot data immediately
	! if (MODULO(i-1,10) == 0) THEN
	! 	 call plot_data_immediately(sstate,i,ptr_gctrl)
	!  end if

	! Write data to file
	call write_data_to_file(sstate,i)

end do


!print out elapsed time:
call system_clock ( t2, clock_rate, clock_max )
write ( *, * ) 'Elapsed time = ', real (t2-t1)/real(clock_rate), 'seconds, which is an average of ' &
														 & ,  real (t2-t1)/real(clock_rate)/real(nFrames), ' seconds per frame.'


! Cleanup
call free_state(sstate)

end program sph

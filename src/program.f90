program sph

	! Main program of project
	
	!Created by Jannik Zuern on 05/16/2016
	!Last modified: 05/16/2016

	use gnufor2       !Module: gnuplot fortran visualizations
	use integrate     !Module: definition of integration functions
	use sphfunctions  !Module: implementation of sph mechanics
	use linkedlists   !Module: Linked lists definition
	use util          !Module: utiliary functions

	implicit none

! 	integer, parameter	:: N1=50
! 	real(kind=8)		:: x1(N1), f1(N1)
! 	integer :: i
! ! generate data for 2D plots
! 	do  i=1,N1
! 		x1(i)=5.0*i/N1
! 	end do
! 	f1=sin(2*x1)
!
!
! 	print *,'Example 1: simple 2D graph'
! 	print *,'call plot(x1,f1)'
! 	call plot(x1,f1)
! 	print *,'press ENTER to go to the next example'
! 	read  *



! ACTUAL PROGRAM STARTING
!-------------------------------------------------------------
!-------------------------------------------------------------

! Define variables
DOUBLE PRECISION, dimension (9) :: parameter


! Parse input parameter file

call parse_input(parameter)

! Initialize sim param vector

! Initialize sim state Class instantiation (= Object)

! ....


! Simulation loop




! Cleanup




end program sph

module util

  ! Auxiliary utility types and subroutines and functions

  !Created by Jannik Zuern on 05/16/2016
	!Last modified: 05/16/2016

implicit none

  type :: systemstate

    ! The systemstate data type hold all arrays of the physical simulation variables
    ! like particle position, velocity, acceleration, density,..

     double precision,allocatable,dimension(:) :: x    ! position
     double precision,allocatable,dimension(:) :: v    ! velocity
     double precision,allocatable,dimension(:) :: vh   ! save temp velocty in leapfrog step
     double precision,allocatable,dimension(:) :: a    ! acceleration
     double precision,allocatable,dimension(:) :: rho  ! density
     double precision,allocatable,dimension(:) :: mass ! mass


  end type systemstate

contains


  subroutine parse_input(parameter)

    ! parse_input writes parameter from parameter input file into parameter array

    DOUBLE PRECISION, DIMENSION(9)  :: parameter
    CHARACTER(len=32) :: filename
    CHARACTER(len=32) :: line
    integer :: i = 1, j


    CALL get_command_argument(i, filename)

    ! process file filename
    print *, "Parsing file ", filename
    open (unit = 100, file = filename, action = 'read')

    do j = 1,  9
      read(100,*) line
      read(line,*) parameter(j)   ! convert from character type to dp
    end do

    close(unit=100)


  end subroutine







end module util

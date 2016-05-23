module util

  ! Auxiliary utility types and subroutines and functions

  !Created by Jannik Zuern on 05/16/2016
	!Last modified: 05/19/2016


implicit none


private
public :: systemstate, alloc_state, parse_input, circ_indicator, box_indicator,free_state

  type systemstate

    ! The systemstate data type hold all arrays of the physical simulation variables
    ! like particle position, velocity, acceleration, density,..

     double precision,allocatable,dimension(:) :: x    ! position
     double precision,allocatable,dimension(:) :: v    ! velocity
     double precision,allocatable,dimension(:) :: vh   ! save temp velocty in leapfrog step
     double precision,allocatable,dimension(:) :: a    ! acceleration
     double precision,allocatable,dimension(:) :: rho  ! density

     double precision                          :: mass ! mass
     integer                                   :: nParticles ! number of particles

  end type systemstate





contains


  function box_indicator (x,y) result(res)
    implicit none
    double precision                :: x,y
    logical                         :: tmp
    integer                         :: res

    res = 0
    ! print *, x , y
    tmp = (x > 0.3) .AND. (y > 0.3) .AND. (x < 0.35) .AND. (y < 0.35)
    if (tmp .eqv. .true.) THEN

      res = 1
    end if

  end function


  integer function circ_indicator (x,y) result(res)
    implicit none
    double precision, intent(in)    :: x,y
    logical                         :: tmp
    ! integer                         :: res
    double precision                :: dx,dy,r2

    dx = x-0.5
    dy = y-0.6
    r2 = dx*dx + dy*dy
    res = 0

    tmp = (r2 > 0.1*0.1) .AND. (r2 < 0.3*0.3)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function



  subroutine alloc_state(state)
    type(systemstate) :: state !system state object

    ! allocate all arrays
    allocate(state%x   (2*state%nParticles))
    allocate(state%v   (2*state%nParticles))
    allocate(state%vh  (2*state%nParticles))
    allocate(state%a   (2*state%nParticles))
    allocate(state%rho (  state%nParticles))
    
    state%x = 0
    state%v = 0
    state%vh = 0
    state%a = 0
    state%rho = 0

  end subroutine


  subroutine free_state(state)
    type(systemstate) :: state !system state object

    ! deallocate all arrays
    deallocate(state%x  )
    deallocate(state%v  )
    deallocate(state%vh )
    deallocate(state%a  )
    deallocate(state%rho)


  end subroutine

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
    ! iterate through each line
    do j = 1,  9
      read(100,*) line
      print *, line
      read(line,*) parameter(j)   ! convert from character type to dp
    end do

    close(unit=100)

    print *, "Parsing completed"
  end subroutine






  subroutine plot_points()

  end subroutine











end module util

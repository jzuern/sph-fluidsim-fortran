module util

  !Auxiliary utility types, subroutines, and functions

  !Created by Jannik Zuern on 05/16/2016
	!Last modified: 05/24/2016


implicit none


private
public :: systemstate, sim_parameter, alloc_state, parse_input, circ_indicator, box_indicator , &
free_state , plot_data_immediately, plot_data_from_file,write_data_to_file , initialize_parameters, Pi

  double precision, PARAMETER :: Pi = 3.1415927d0

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

  type sim_parameter

    ! The sim_parameter data type holds all parameter that are relevant for
    ! the simulation. I.e. particle size, time step size, number of frames to
    ! calculate, etc.


     integer          :: nframes            ! Number of frames
     integer          :: nSteps_per_frame   ! Number of steps per frame. Default: 100
     double precision :: h                  ! Size of particles (radius). Also: Distance of particles in initial configuration
     double precision :: dt                 ! Time step size Default: 1E-4
     double precision :: rho0               ! Reference density Default: 1000
     integer          :: k                  ! Bulk modulus (Kompressionsmodul) Default: 1E3
     double precision :: mu                 ! Viscosity Default: 0.1
     double precision :: g                  ! gravity strength Default: 9.81
     double precision :: rcut               ! Cutoff radius (in multiples of total size of simulation grid for linked lists neighbor tracking

  end type sim_parameter



contains


  function box_indicator (x,y) result(res)
    implicit none
    double precision                :: x,y
    logical                         :: tmp
    integer                         :: res

    res = 0
    ! print *, x , y
    tmp = (x > 0.3d0) .AND. (y > 0.0d0) .AND. (x < 0.0d0) .AND. (y < 0.7d0)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function


  function circ_indicator (x,y) result(res)
    implicit none
    double precision, intent(in)    :: x,y
    logical                         :: tmp
    integer                         :: res
    double precision                :: dx,dy,r2

    dx = x-0.2d0
    dy = y-0.2d0
    r2 = dx*dx + dy*dy
    res = 0

    tmp = (r2 > 0.1d0*0.1d0) .AND. (r2 < 0.5d0*0.5d0)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function



  subroutine alloc_state(state,params)
    type(systemstate) :: state !system state object
    type(sim_parameter)											:: params
    double precision :: rho0

    ! allocate all arrays
    allocate(state%x   (2*state%nParticles))
    allocate(state%v   (2*state%nParticles))
    allocate(state%vh  (2*state%nParticles))
    allocate(state%a   (2*state%nParticles))
    allocate(state%rho (  state%nParticles))

    state%x       = 0.d0
    state%v       = 0.d0
    state%vh      = 0.d0
    state%a       = 0.d0
    state%rho     = params%rho0
    state%mass    = 1.0d0

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


  subroutine write_data_to_file(sstate,i)

    type(systemstate)                           :: sstate
    integer                                     :: n,i,k
    double precision, allocatable, dimension(:) :: x,y
    character(len=100) :: filename
    character(len=5)   :: dummy
    character(len=8)   :: fmt ! format descriptor

    fmt = '(I5.5)'   ! an integer of width 5 with zeros at the left

    write (dummy,fmt) i ! converting integer to string using a 'internal file'
    filename='data/frame'//trim(dummy)//'.dat'
    open(unit = 10, file = filename)

    n = sstate%nParticles
    x = sstate%x(1:2*n:2)
    y = sstate%x(2:2*n:2)

    do k = 1,n
      write (10,*) x(k),y(k)
    end do


    close(10)
  end subroutine


  subroutine plot_data_from_file(sstate,i)
    use gnufor2
    type(systemstate)                           :: sstate
    integer                                     :: i

      character(len=100) :: datafile
      character(len=100) :: commandfile
      character(len=5)   :: dummy
      character(len=8)   :: fmt ! format descriptor
      character(len=100) :: line1,line2


      fmt = '(I5.5)'   ! an integer of width 5 with zeros at the left
      write (dummy,fmt) i ! converting integer to string using a 'internal file'
      datafile = "dataframe" // trim(dummy) // ".dat"


      !writing command file (here we specify the command file, from where gnuplot will
      !                      get all the info about plotting specifications)
      commandfile = "command_file.txt"

      open (unit=11, file=commandfile, status="replace")
      write ( 11,*) "set terminal  x11 size 600,600"
      line1 = "set output '" // datafile // "'" !TODO: find out why linebreak in string
      write ( 11,*) line1
      write ( 11,*) "unset key"
      write ( 11,*) "set xrange [0:1]"
      write ( 11,*) "set yrange [0:1]"
      write ( 11,*) "set grid"
      line2 = 'plot "' // commandfile // '" using 1:2 with points pointtype  4 linecolor rgb "blue" linewidth 1'
      write ( 11,*) line2 !TODO: find out why linebreak in string
      write ( 11,*) "pause 0.100E+00"
      write ( 11,*) "q"

      !calling gnuplot
      call run_gnuplot (commandfile)

  end subroutine


  subroutine plot_data_immediately(sstate,i)
    use gnufor2

    type(systemstate)                           :: sstate
    integer                                     :: n,i
    double precision, allocatable, dimension(:) :: x,y

    character(len=100) :: file
    character(len=5)   :: dummy
    character(len=8)   :: fmt ! format descriptor

    n = sstate%nParticles
    x = sstate%x(1:2*n:2)
    y = sstate%x(2:2*n:2)

    fmt = '(I5.5)'   ! an integer of width 5 with zeros at the left
    write (dummy,fmt) i ! converting integer to string using a 'internal file'
    file='data/frame'//trim(dummy)//'.dat'

    call plot(x,y," 4.",pause=0.1, persist = "yes", terminal=" x11 size 600,600", filename = file)

  end subroutine



  function initialize_parameters (params) result(param_type)
    implicit none
    DOUBLE PRECISION, DIMENSION(9)             :: params
    type(sim_parameter)                        :: param_type

    param_type%nframes           = params(1)
    param_type%nSteps_per_frame  = params(2)
    param_type%h                 = params(3)
    param_type%dt                = params(4)
    param_type%rho0              = params(5)
    param_type%k                 = params(6)
    param_type%mu                = params(7)
    param_type%g                 = params(8)
    param_type%rcut              = params(9)


  end function






end module util

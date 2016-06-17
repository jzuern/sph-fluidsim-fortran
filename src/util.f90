module util

  !Auxiliary utility types, subroutines, and functions

  !Created by Jannik Zuern on 05/16/2016
	!Last modified: 06/08/2016


implicit none


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
     integer                                   :: nLiquidParticles ! number of particles in liquid
     integer                                   :: nSolidParticles  ! number of particles in solid object
     integer                                   :: nParticles       ! total number of particles in simulation

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
     double precision :: rcut_x             ! Cutoff radius in x-direction for cells (as fraction of total size of simulation grid)
     double precision :: rcut_y             ! Cutoff radius in y-direction for cells (as fraction of total size of simulation grid)
     double precision :: dphi               ! rotation speed of cross

  end type sim_parameter



contains


  function cross_indicator (x,y) result(res)

    !indicates whether point (x,y) lies within a cross (res == 1) or not (res == 0)

    implicit none

    double precision,intent(in)     :: x,y
    logical                         :: tmp,b1,b2
    integer                         :: res

    double precision :: diameter  = 0.7d0  ! diameter of rotating cross
    double precision :: thickness = 0.02d0 ! thickness of roating cross beams
    double precision, dimension(2) :: center = (/0.5d0, 0.5d0/) ! rotation center


    res = 0

    ! check wheter x and y coordinates lay within domain of cross
    b1 = ( ABS(x - center(1)) < diameter/2  ) .AND. ( ABS(y - center(2)) < thickness/2  )
    b2 = ( ABS(y - center(2)) < diameter/2  ) .AND. ( ABS(x - center(1)) < thickness/2  )

    if (b1 .OR. b2) THEN
      res = 1
    end if

  end function





  function box_indicator (x,y) result(res)

    !indicates whether point (x,y) lies within rectangular box (res == 1) or not (res == 0)

    implicit none

    double precision,intent(in)     :: x,y
    logical                         :: tmp
    integer                         :: res

    res = 0
    tmp = (x < 1.0d0) .AND. (x > 0.5d0) .AND. (y > 0.8d0) .AND. (y < 1.0d0)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function


  function circ_indicator (x,y) result(res)

    !indicates whether point (x,y) lies within circle (res == 1) or not (res == 0)

    implicit none

    double precision, intent(in)    :: x,y
    logical                         :: tmp
    integer                         :: res
    double precision                :: dx,dy,r2

    dx = x-0.5d0
    dy = y-0.5d0
    r2 = dx*dx + dy*dy
    res = 0

    tmp = (r2 > 0.2d0*0.2d0) .AND. (r2 < 0.6d0*0.6d0)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function


  subroutine alloc_state(state,params)

    ! allocates memory for system state variables

    type(systemstate)                       :: state !system state object
    type(sim_parameter)											:: params!simulation parameter object
    double precision                        :: rho0

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

    ! deallocates variables and frees memory

    type(systemstate) :: state !system state object

    ! deallocate all arrays
    deallocate(state%x  )
    deallocate(state%v  )
    deallocate(state%vh )
    deallocate(state%a  )
    deallocate(state%rho)

  end subroutine



  subroutine initialize_parameters (params)

    ! reads sim_parameter.dat file using the FORTRAN namelist feature

    implicit none
    type(sim_parameter)                        :: params
    CHARACTER(len=32)                          :: inputfile

    integer          :: nframes
    integer          :: nSteps_per_frame
    double precision :: h
    double precision :: dt
    double precision :: rho0
    double precision :: k
    double precision :: mu
    double precision :: g
    double precision :: rcut_x
    double precision :: rcut_y
    double precision :: dphi


    namelist /SIMPARAMETER/nframes,nSteps_per_frame,h,dt,rho0,k,mu,g,rcut_x,rcut_y,dphi

    CALL get_command_argument(1, inputfile)
    print *, "Parsing file ", inputfile


    open(unit=10, file=inputfile)
    read(10,NML=SIMPARAMETER)
    close(10)

    params%nframes = nframes
    params%nSteps_per_frame = nSteps_per_frame
    params%h = h
    params%dt = dt
    params%rho0 = rho0
    params%k = k
    params%mu = mu
    params%g = g
    params%rcut_x = rcut_x
    params%rcut_y = rcut_y
    params%dphi = dphi

    print *, "...Parsing completed "


  end subroutine


  subroutine write_data_to_file(sstate,i)

    ! writes particle data for each frame to files

    type(systemstate),intent(in)                :: sstate
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

    ! writes particle positions into file
    do k = 1,n
      write (10,*) x(k),y(k)
    end do


    close(10)

  end subroutine


  subroutine plot_data_from_file(sstate,i)

      ! plots particles from particle data files

      use gnufor2
      type(systemstate),intent(in)                :: sstate
      integer                                     :: i

      character(len=100) :: datafile
      character(len=100) :: commandfile
      character(len=5)   :: dummy
      character(len=8)   :: fmt ! format descriptor
      character(len=100) :: line1,line2


      fmt = '(I5.5)'   ! an integer of width 5 with zeros at the left
      write (dummy,fmt) i ! converting integer to string using a 'internal file'
      datafile = "data/frame" // trim(dummy) // ".dat"

      !writing command file (here we specify the command file, from where gnuplot will
      !                      get all the info about plotting specifications)
      commandfile = "command_file.txt"

      open (unit=11, file=commandfile, status="replace")
      write ( 11,'(a,i2,a)') "set terminal  x11 size 1000,1000"
      write ( 11,'(a,i2,a)') "unset key"
      write ( 11,'(a,i2,a)') "set xrange [0:1]"
      write ( 11,'(a,i2,a)') "set yrange [0:1]"
      write ( 11,'(a,i2,a)') "set grid"
      write ( 11,'(a,i2,a)') 'plot "' // trim (datafile) //'" using 1:2 with points pointtype 65 linecolor rgb "blue" linewidth 1'
      write ( 11,'(a,i2,a)') "pause 0.100E+00" ! pause of 0.1 s
      write ( 11,'(a,i2,a)') "q"
      close (11)

      !calling gnuplot
      call run_gnuplot (commandfile)

  end subroutine


  subroutine plot_data_immediately(sstate,i)

    ! plots data without saving to file

    use gnufor2

    type(systemstate),intent(in)                :: sstate
    integer                                     :: n,i
    double precision, allocatable, dimension(:) :: x,y

    character(len=100) :: file
    character(len=5)   :: dummy
    character(len=8)   :: fmt ! format descriptor

    n = sstate%nParticles
    x = sstate%x(1:2*n:2)
    y = sstate%x(2:2*n:2)

    fmt = '(I5.5)'                         ! an integer of width 5 with zeros at the left
    write (dummy,fmt) i                    ! converting integer to string using a 'internal file'
    file='data/frame'//trim(dummy)//'.dat' ! concatenating stuff to create file name

    ! actual plotting command:
    call plot(x,y," 5.",pause=0.5, persist = "yes", terminal=" x11 size 800,800", filename = file)

  end subroutine





end module util

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
     logical          :: mill               ! Decide whether a watermill is in the computational domain or not
     double precision :: dphi               ! rotation speed of watermill


  end type sim_parameter



contains


  function cross_indicator (x,y,params) result(res)

    !indicates whether point (x,y) lies within a cross (res == 1) or not (res == 0)

    implicit none

    type(sim_parameter)             :: params
    double precision,intent(in)     :: x,y
    logical                         :: b1,b2
    integer                         :: res

    double precision :: diameter  = 0.7d0  ! diameter of rotating cross
    double precision :: thickness = 0.02d0 ! thickness of roating cross beams
    double precision, dimension(2) :: center = (/0.5d0, 0.5d0/) ! rotation center



    if(params%mill) then

      res = 0

      ! check wheter x and y coordinates lay within domain of cross
      b1 = ( ABS(x - center(1)) < diameter/2  ) .AND. ( ABS(y - center(2)) < thickness/2  )
      b2 = ( ABS(y - center(2)) < diameter/2  ) .AND. ( ABS(x - center(1)) < thickness/2  )

      if (b1 .OR. b2) THEN
        res = 1
      end if

    else
      res = 0
    end if

  end function



  function box_indicator (x,y) result(res)

    !indicates whether point (x,y) lies within rectangular box (res == 1) or not (res == 0)

    implicit none

    double precision,intent(in)     :: x,y
    logical                         :: tmp
    integer                         :: res
    double precision                :: xmin,xmax,ymin,ymax

    xmin = 0.5d0
    xmax = 1.0d0
    ymin = 0.8d0
    ymax = 1.0d0

    res = 0
    tmp = (x < xmax) .AND. (x > xmin) .AND. (y > ymin) .AND. (y < ymax)
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
    double precision                :: dx,dy,r2,x_offset,y_offset,rmin,rmax

    x_offset = 0.3d0 ! x-coordinates of blob center
    y_offset = 0.3d0 ! y-coordinates of blob center

    rmin = 0.0d0 ! inner radius of circular blob
    rmax = 0.2d0 ! outer radius of circular blob

    dx = x - x_offset
    dy = y - y_offset
    r2 = dx*dx + dy*dy
    res = 0

    tmp = (r2 > rmin*rmin) .AND. (r2 < rmax*rmax)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function


  subroutine alloc_state(state,params)

    ! allocates memory for system state variables

    type(systemstate)                       :: state !system state object
    type(sim_parameter)											:: params!simulation parameter object

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
    logical          :: mill
    double precision :: dphi


    namelist /SIMPARAMETER/nframes,nSteps_per_frame,h,dt,rho0,k,mu,g,rcut_x,rcut_y,mill,dphi

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
    params%mill   = mill
    params%dphi = dphi

    print *, "...Parsing completed "

  end subroutine




  subroutine write_data_to_file(sstate,i)

    ! writes particle data for each frame to files

    type(systemstate),intent(in)                :: sstate
    integer                                     :: n,i,k
    double precision, allocatable, dimension(:) :: x,y,rho
    character(len=100) :: filename
    character(len=5)   :: dummy
    character(len=8)   :: fmt ! format descriptor

    ! fmt = '(I5.5)'   ! an integer of width 5 with zeros at the left
    fmt = '(I0)'   ! an integer of width 5 with zeros at the left


    write (dummy,fmt) i ! converting integer to string using a 'internal file'
    ! write(dummy,*) i
    filename='data/frame'//trim(dummy)//'.dat'


    open(unit = 10, file = filename)

    n = sstate%nParticles
    x = sstate%x(1:2*n:2)
    y = sstate%x(2:2*n:2)
    rho = sstate%rho(1:n:1)


    ! writes particle positions into file
    do k = 1,n
      write (10,*) x(k),y(k),rho(k)
    end do


    close(10)

  end subroutine


  subroutine plot_data_immediately(sstate,i,ptr_gctrl)

    ! plots data without saving to file (for debugging / quick glance at results)

    ! Attention: only works for the first 64 frames since the gnuplotfortran library
    !            for plotting only allows for 64 temp files to exist at the same time
    !            (this is due to a known bug, which has never been fixed by the authors..)


    use datatypes, only : i4b
    use gnuplot_module_data, only : gnuplot_ctrl
    use gnuplot_module

    type(systemstate),intent(in)                :: sstate
    integer :: i
    type(gnuplot_ctrl), pointer :: ptr_gctrl
    integer(i4b) :: status
    integer(i4b) :: n

    n = sstate%nParticles

    status = gnuplot_resetsession(ptr_gctrl)
    ! print *, status
    status = gnuplot_plot2d(ptr_gctrl,n,sstate%x(1:2*n:2),sstate%x(2:2*n:2),'particles')

  end subroutine





  subroutine invoke_gnuplot(ptr_gctrl)

            use datatypes, only : i4b,dp,lgc
            use gnuplot_module_data, only : PI_D,gnuplot_ctrl,GNUPLOT_SHOWDEBUG,GNUPLOT_SHOWWARNINGS,GP_CMD_SIZE
            use gnuplot_module

            integer(i4b), external :: fortran_getchar

            integer(i4b) :: numpoints=50,ii=0
            integer(i4b) :: status=0
            type(gnuplot_ctrl), pointer :: ptr_gctrl

            character(len=1) :: debug

            ! disable Gnuoplot warnings
            GNUPLOT_SHOWWARNINGS=.false.

            ! invokng gnuplot session
            ptr_gctrl=>gnuplot_init('-persist')
            if(.not.associated(ptr_gctrl)) stop 'Failed : to initiate a gnuplot session'

            print*, 'Press Enter to continue ...'

            status=fortran_getchar(debug)

            status=gnuplot_setrange(ptr_gctrl,'x',0.0_dp,1.0_dp)
            status=gnuplot_setrange(ptr_gctrl,'y',0.0_dp,1.0_dp)
            status=gnuplot_settitle(ptr_gctrl,'Particle visualization')
            status=gnuplot_setaxislabel(ptr_gctrl,'x','x')
            status=gnuplot_setaxislabel(ptr_gctrl,'y','y')
            status=gnuplot_setscale(ptr_gctrl,'x','NLG')
            status=gnuplot_setscale(ptr_gctrl,'y','NLG')
            status=gnuplot_setstyle(ptr_gctrl,'points')

      end subroutine invoke_gnuplot



end module util

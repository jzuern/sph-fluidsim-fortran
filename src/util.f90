module util

  !Auxiliary utility types, subroutines, and functions

  !Created by Jannik Zuern on 05/16/2016
	!Last modified: 07/05/2016


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
     double precision                          :: energy ! potential energy and energy of motion
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
     double precision :: rcut_z             ! Cutoff radius in z-direction for cells (as fraction of total size of simulation grid)
     logical          :: mill               ! Decide whether a watermill is in the computational domain or not
     double precision :: dphi               ! rotation speed of watermill


  end type sim_parameter



contains


  subroutine update_solid_particles_positions(sstate,params)

    ! updates solid particles position solely based on passed time

    type(systemstate)                           :: sstate !system state object
    type(sim_parameter)											    :: params
    double precision, dimension(3)              :: center = (/0.5d0, 0.5d0 , 0.5d0/) ! rotation center
    integer                                     :: particle
    double precision                            :: phi,phi_new, radius,dist_x,dist_y,dist_z

    do particle = 1 , sstate%nSolidParticles

      ! convert positions from cartesian to polar coordinates on order to make easy calculations of the new particle positions
      dist_x = sstate%x(3*particle-2) - center(1)
      dist_z = sstate%x(3*particle-0) - center(2)

      radius  =  SQRT (dist_x*dist_x + dist_z*dist_z )
      phi     =  ATAN2 (dist_x , dist_z)

      phi_new =  phi + params%dphi


      ! convert positions back to cartesian coordinates
      sstate%x(3*particle-2) = center(1) + radius * SIN(phi_new)
      sstate%x(3*particle-0) = center(2) + radius * COS(phi_new)
      ! sstate%x(3*particle-0) = sstate%x(3*particle-0) ! y-coordinate stays the same



    end do

  end subroutine



  function cross_indicator (x,y,z,params) result(res)

    !indicates whether point (x,y) lies within a cross (res == 1) or not (res == 0)

    implicit none

    type(sim_parameter)             :: params
    double precision,intent(in)     :: x,y,z
    logical                         :: b1,b2,b3
    integer                         :: res

    double precision :: diameter  = 0.70d0  ! diameter of rotating cross
    double precision :: thickness = 0.05d0 ! thickness of roating cross beams
    double precision, dimension(2) :: center = (/0.5d0, 0.5d0/) ! rotation center (x.y coordinates)



    if(params%mill) then

      res = 0

      ! check whether x and z coordinates lay within domain of cross
      b1 = ( ABS(x - center(1)) < diameter/2  ) .AND. ( ABS(z - center(2)) < thickness/2  )
      b2 = ( ABS(z - center(2)) < diameter/2  ) .AND. ( ABS(x - center(1)) < thickness/2  )

      ! set cross height
      b3 = y < 0.5

      if ((b1 .OR. b2) .AND. b3) THEN
        res = 1
      end if

    else
      res = 0
    end if

  end function



  function box_indicator (x,y,z) result(res)

    !indicates whether point (x,y) lies within rectangular box (res == 1) or not (res == 0)

    implicit none

    double precision,intent(in)     :: x,y,z
    logical                         :: tmp
    integer                         :: res
    double precision                :: xmin,xmax,ymin,ymax,zmin,zmax

    xmin = 0.5d0
    xmax = 1.0d0
    ymin = 0.8d0
    ymax = 1.0d0
    zmin = 0.0d0
    zmax = 0.2d0

    res = 0
    tmp = (x < xmax) .AND. (x > xmin) .AND. (y > ymin) .AND. (y < ymax) .AND. (z > zmin) .AND. (z < zmax)
    if (tmp .eqv. .true.) THEN
      res = 1
    end if

  end function


  function circ_indicator (x,y,z) result(res)

    !indicates whether point (x,y) lies within circle (res == 1) or not (res == 0)

    implicit none

    double precision, intent(in)    :: x,y,z
    logical                         :: tmp
    integer                         :: res
    double precision                :: dx,dy,dz,r2,x_offset,y_offset,z_offset,rmin,rmax

    x_offset = 0.8d0 ! x-coordinates of blob center
    y_offset = 0.3d0 ! y-coordinates of blob center
    z_offset = 0.5d0 ! z-coordinates of blob center

    rmin = 0.00d0 ! inner radius of circular blob
    rmax = 0.20d0 ! outer radius of circular blob

    dx = x - x_offset
    dy = y - y_offset
    dz = z - z_offset
    r2 = dx*dx + dy*dy + dz*dz
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
    allocate(state%x   (3*state%nParticles))
    allocate(state%v   (3*state%nParticles))
    allocate(state%vh  (3*state%nParticles))
    allocate(state%a   (3*state%nParticles))
    allocate(state%rho (1*state%nParticles))

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
    double precision :: rcut_z
    logical          :: mill
    double precision :: dphi


    namelist /SIMPARAMETER/nframes,nSteps_per_frame,h,dt,rho0,k,mu,g,rcut_x,rcut_y,rcut_z,mill,dphi

    CALL get_command_argument(1, inputfile)
    print *, "Parsing file ", inputfile


    open(unit=10, file=inputfile)
    read(10,NML=SIMPARAMETER)
    close(10)

    params%nframes              = nframes
    params%nSteps_per_frame     = nSteps_per_frame
    params%h                    = h
    params%dt                   = dt
    params%rho0                 = rho0
    params%k                    = k
    params%mu                   = mu
    params%g                    = g
    params%rcut_x               = rcut_x
    params%rcut_y               = rcut_y
    params%rcut_z               = rcut_z
    params%mill                 = mill
    params%dphi                 = dphi

    print *, "...Parsing completed "

  end subroutine




  subroutine write_data_to_file(sstate,i)

    ! writes particle data for each frame to files

    type(systemstate),intent(in)                :: sstate
    integer                                     :: n,i,k
    double precision, allocatable, dimension(:) :: x,y,z,rho
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
    x = sstate%x(1:3*n:3)
    y = sstate%x(2:3*n:3)
    z = sstate%x(3:3*n:3)
    rho = sstate%rho(1:n:1)


    ! writes particle positions into file
    do k = 1,n
      write (10,*) x(k),z(k),y(k),rho(k)
    end do


    close(10)

  end subroutine


  subroutine plot_data_immediately(sstate,i,ptr_gctrl)

    ! plots data without saving to file (for debugging / quick glance at results)

    ! Attention: only works for the first 64 frames since the gnuplotfortran library
    !            for plotting only allows for 64 temp files to exist at the same time
    !            (this is due to a known bug, which has never been fixed by the  library authors...)


    use datatypes, only : i4b
    use gnuplot_module_data, only : gnuplot_ctrl
    use gnuplot_module


    type(systemstate),intent(in)                   :: sstate
    integer                                        :: i,iter,c
    type(gnuplot_ctrl), pointer                    :: ptr_gctrl
    integer(i4b)                                   :: status
    integer(i4b)                                   :: n,n_c
    double precision, allocatable, dimension (:,:) :: y
    double precision, allocatable, dimension(:)    :: x_coords,z_coords

    n = sstate%nParticles

    c = 5 ! plot only every c-th particle (it is more efficient)
    n_c = (n/c)

    allocate (y(n_c,n_c))
    allocate (x_coords(n_c))
    allocate (z_coords(n_c))

    ! preparing data structure for plotting
    x_coords = sstate%x(1:3*n:3*c)
    z_coords = sstate%x(3:3*n:3*c)

    y = -1.d0 ! initialize z to value of -1.0 so that corresponding particles won't be visible in plot
    do iter = 1,n_c
      y(iter,iter) = sstate%x(3*iter*c - 1)
    end do


    status = gnuplot_resetsession(ptr_gctrl) ! remove previous plots from current Gnuplot session
    if(status.ne.0) stop 'Failed to reset Gnuplot session.'

    status = gnuplot_plot3d(ptr_gctrl , n_c , n_c , x_coords , z_coords , y)
    if(status.ne.0) stop 'Failed : to splot (3D) in subroutine plot_data_immediately'


  end subroutine





  subroutine invoke_gnuplot(ptr_gctrl)

    ! invoke a Gnuplot session from within FORTRAN

      use datatypes, only : i4b,dp,lgc
      use gnuplot_module_data, only : gnuplot_ctrl,GNUPLOT_SHOWWARNINGS
      use gnuplot_module

      integer(i4b), external          :: fortran_getchar
      integer(i4b)                    :: status=0
      type(gnuplot_ctrl), pointer     :: ptr_gctrl
      character(len=1)                :: debug

      ! disable Gnuoplot warnings
      GNUPLOT_SHOWWARNINGS=.false.

      ! invokng gnuplot session
      ptr_gctrl=>gnuplot_init('-persist')
      if(.not.associated(ptr_gctrl)) stop 'Failed to initiate a gnuplot session.'

      print*, '  '
      print*, 'Press Enter to continue ...'

      status=fortran_getchar(debug)

      status=gnuplot_set(ptr_gctrl,'terminal x11 size 800,800')

      status=gnuplot_setrange(ptr_gctrl,'x',0.0_dp,1.0_dp)
      status=gnuplot_setrange(ptr_gctrl,'y',0.0_dp,1.0_dp)
      status=gnuplot_setrange(ptr_gctrl,'z',0.0_dp,1.0_dp)
      status=gnuplot_settitle(ptr_gctrl,'Particle visualization')
      status=gnuplot_setaxislabel(ptr_gctrl,'x','x')
      status=gnuplot_setaxislabel(ptr_gctrl,'y','y')
      status=gnuplot_setaxislabel(ptr_gctrl,'z','z')
      ! status=gnuplot_setscale(ptr_gctrl,'x','NLG')
      ! status=gnuplot_setscale(ptr_gctrl,'y','NLG')
      ! status=gnuplot_setscale(ptr_gctrl,'z','NLG')
      status=gnuplot_setstyle(ptr_gctrl,'points')



  end subroutine invoke_gnuplot



end module util

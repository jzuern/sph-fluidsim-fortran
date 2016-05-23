module sphfunctions

  ! Implementation of functions and subroutines for the SPH method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/19/2016


implicit none


private

public :: reflect_bc, compute_density_with_ll, compute_accel,check_state, init_particles

contains


  subroutine init_particles(sstate,params)

    use util
    type(systemstate)                  :: sstate !system state object
    double precision, dimension(9)     :: params

    call place_particles(sstate,params)
    ! call normalize_mass(sstate,params)

    print *, "init particles completed"
  end subroutine



  subroutine check_state(sstate)
    use util
    type(systemstate) :: sstate !system state object

    integer :: i
    double precision :: xi, yi

    do i = 1,sstate%nParticles
      xi = sstate%x(2*i - 1)
      yi = sstate%x(2*i - 0)
      ! assert statements unavailable in F90
      !         assert( xi >= 0 || xi <= 1 );
      !         assert( yi >= 0 || yi <= 1 );
    end do



  end subroutine




  subroutine reflect_bc(sstate)
    use util
    type (systemstate) :: sstate




  end subroutine



  subroutine normalize_mass(sstate,params)
    use util
    type (systemstate)             :: sstate
      DOUBLE PRECISION, DIMENSION(9)  :: params
    double precision               :: rho0   ! reference density
    double precision               :: rho2s,rhos
    integer                        :: i

    rhos = 0.d0
    rho2s = 0.d0
    rho0 = params(5)

    sstate%mass = 1

    do i = 1,sstate%nParticles
      rho2s = rho2s + sstate%rho(i)*sstate%rho(i)
      rhos  = rhos + sstate%rho(i)
    end do

    sstate%mass = sstate%mass * (rho0*rhos / rho2s)



  end subroutine



  subroutine compute_density_with_ll(sstate, params, ll, lc)
    use util
    type (systemstate)             :: sstate
    DOUBLE PRECISION, DIMENSION(9)  :: params
    integer, dimension(:)          :: ll
    integer, dimension(:,:)        :: lc


    ! TODO: implement body




  end subroutine


  subroutine compute_accel(s,params,ll,lc)
    use util
    type (systemstate)           :: s
    DOUBLE PRECISION, DIMENSION(9)  :: params
    integer, dimension(:)        :: ll
    integer, dimension(:,:)      :: lc

    ! TODO: implement body

  end subroutine


  subroutine place_particles(sstate,params)

    use util
    type (systemstate)                    :: sstate
    DOUBLE PRECISION, DIMENSION(9)        :: params
    double precision                      :: h,hh
    integer                               :: count, p
    double precision                      :: x,y, rd

    h = params(3) ! size of particles
    hh = h/1.0d0  ! why?


    count = 0


    x = 0.d0
    do while (x < 1.d0)
      y = 0.d0
      do while (y < 1.d0)
        count = count + box_indicator(x,y)
        y = y + hh
      end do
      x = x + hh
    end do

    print *, "Number of particles in simulation: " , count

    sstate%nParticles = count     ! set number of particles to counted value
    call alloc_state(sstate)      ! allocate fields for x,v,vh,a

    p = 1 ! particle iterator

    x = 0.d0
    do while (x < 1.d0)
      y = 0.d0
      do while (y < 1.d0)

        if (box_indicator(x,y) /= 0) THEN
          ! CALL init_random_seed()         ! TODO: do I need initialization of random seed?
          CALL RANDOM_NUMBER(rd)            ! random number between 0 and 1
          rd = rd * 0.000001d0              ! TODO: what is correct size for this?
          print *, x , y , rd
          sstate%x(2*p-1) = x;
          sstate%x(2*p-0) = y;
          sstate%v(2*p-1) = rd; ! initialize with small initial velocity (-> no symmetries in behaviour)
          sstate%v(2*p-0) = rd;
          p = p + 1
        end if

        y = y + hh
      end do
      x = x + hh
    end do

    print *, "placed particles"


  end subroutine





end module sphfunctions

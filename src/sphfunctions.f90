module sphfunctions

  ! Implementation of functions and subroutines for the SPH method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/19/2016


implicit none


private

public :: reflect_bc, damp_reflect, compute_density_with_ll, compute_accel,check_state, init_particles

contains


  subroutine init_particles(sstate,params)

    use util
    type(systemstate)                  :: sstate !system state object
    double precision, dimension(9)     :: params

    call place_particles(sstate,params)
    call normalize_mass(sstate,params)

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
    type (systemstate)                          :: sstate
    DOUBLE PRECISION, DIMENSION(9)              :: params
    integer, allocatable, dimension(:)          :: ll
    integer, allocatable, dimension(:,:)        :: lc
    integer, dimension(4)                       :: ndx,ndy
    integer, dimension(2)                       :: nmax

    integer              :: n, nx ,ny
    integer              :: n1,n2,no
    integer              :: i,j
    double precision     :: mass,h,h2,h8,C,rcut

    real :: dx,dy,r2,z, rho_ij

    h = params(3)
    h2 = h*h
    h8 = h2*h2*h2*h2

    n          = sstate%nParticles
    mass       = sstate%mass

    ndx = (/1,1,0,-1 /)
    ndy = (/0,1,1, 1 /)

    rcut = params(9)             ! is 9th element in sim_param vector....
    nmax(1) = int(floor(1/rcut)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1/rcut)) ! maximum number of cells in each dimension


    do i = 1,nmax(1)
      do j = 1,nmax(2)

        if (lc(i,j) /= -1) THEN
          n1 = lc(i,j)

          do while (n1 /= -1)
            n2 = ll(n1)
            sstate%rho(n1) = sstate%rho(n1) + 4*mass/Pi/h2

            do while(n2 /= -1)
              dx = sstate%x(2*n1-1) - sstate%x(2*n2-1)
              dy = sstate%x(2*n1-0) - sstate%x(2*n2-0)
              r2 = dx*dx + dy*dy
              z  = h2 - r2

              if (z > 0.d0) then
                rho_ij = C*z*z*z
                sstate%rho(n1) = sstate%rho(n1) + rho_ij
                sstate%rho(n2) = sstate%rho(n2) + rho_ij
              end if

              n2 = ll(n2)
            end do

            ! Now the neighboring cells of cell i,j
            do no = 1,4
              nx = i+ndx(no)
              ny = j+ndy(no)
              !boundary conditions
              if (nx <         0) continue
              if (nx > nmax(1)-1) continue
              if (ny <         0) continue
              if (ny > nmax(2)-1) continue

              n2 = lc(nx,ny)

              do while (n2 /= -1)
                dx = sstate%x(2*n1-1) - sstate%x(2*n2-1)
                dy = sstate%x(2*n1-0) - sstate%x(2*n2-0)
                r2 = dx*dx + dy*dy
                z  = h2 - r2
                if (z > 0.d0) then
                  rho_ij = C*z*z*z
                  sstate%rho(n1) = sstate%rho(n1) + rho_ij
                  sstate%rho(n2) = sstate%rho(n2) + rho_ij
                end if
                n2 = ll(n2)
              end do
            end do
            n1 = ll(n1)
          end do
        end if
      end do
    end do

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


  subroutine damp_reflect()


  end subroutine



    subroutine reflect_bc(sstate)
      use util
      type (systemstate) :: sstate

      ! void reflect_bc(sim_state_t* s){
      !
      !         float*  vh = s->vh;
      !         float*  v = s->v;
      !         float*  x = s->x;
      !         int n = s->n;
      !         for (int i = 0; i < n; ++i, x += 2, v += 2, vh += 2) {
      !                 if (x[0] < XMIN) damp_reflect(0, XMIN, x, v, vh);
      !                 if (x[0] > XMAX) damp_reflect(0, XMAX, x, v, vh);
      !                 if (x[1] < YMIN) damp_reflect(1, YMIN, x, v, vh);
      !                 if (x[1] > YMAX) damp_reflect(1, YMAX, x, v, vh);
      !         }
      ! }
      integer :: i,n
      ! Boundaries of computational domain
      double precision :: xmin = 0.0
      double precision :: xmax = 1.0
      double precision :: ymin = 0.0
      double precision :: ymax = 1.0

      n = sstate%nParticles
      do i = 1,n
        if (sstate%x(1) < xmin) then
          call damp_reflect(0,xmin,x,v,vh)
        end if
        if (sstate%x(1) > xmax) then
          call damp_reflect()
        end if
        if (sstate%x(2) < ymin) then
          call damp_reflect()
        end if
        if (sstate%x(2) > ymax) then
          call damp_reflect()
        end if
      end do

    end subroutine




  subroutine compute_accel(s,params,ll,lc)
    use util
    type (systemstate)           :: s
    DOUBLE PRECISION, DIMENSION(9)  :: params
    integer, dimension(:)        :: ll
    integer, dimension(:,:)      :: lc

    ! TODO: implement body




  end subroutine





end module sphfunctions

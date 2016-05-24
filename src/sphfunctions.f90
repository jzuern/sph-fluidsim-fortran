module sphfunctions

  ! Implementation of functions and subroutines for the SPH method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/24/2016


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

    sstate%mass = 1.d0

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

    double precision :: dx,dy,r2,z, rho_ij

    h = params(3)
    h2 = h*h
    h8 = h2*h2*h2*h2

    n          = sstate%nParticles
    mass       = sstate%mass
    !
    ! print *, "in compute_density_with_ll"
    ! print *, sstate%rho
    ! print *, " h2 = " , h2, mass, Pi

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
    hh = h/1.1d0  ! why?


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
    call alloc_state(sstate,params)      ! allocate fields for x,v,vh,a

    p = 1 ! particle iterator

    x = 0.d0
    do while (x < 1.d0)
      y = 0.d0
      do while (y < 1.d0)

        if (box_indicator(x,y) /= 0) THEN
          ! CALL init_random_seed()         ! TODO: do I need initialization of random seed?
          CALL RANDOM_NUMBER(rd)            ! random number between 0 and 1
          rd = rd * 0.01d0              ! TODO: what is correct size for this?

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

  end subroutine


  subroutine damp_reflect(i,which, barrier, sstate)
    use util
    type (systemstate)  :: sstate
    integer             :: which,i
    double precision    :: barrier
    double precision    :: damp, tbounce

    !ignore degenerate cases
    ! if (sstate%v(i+which) == 0) then
    !   return
    ! end if

    ! static void damp_reflect(int which, float barrier, float* x, float* v, float* vh) {

    !         x[0] -= v[0]*(1-DAMP)*tbounce;
    !         x[1] -= v[1]*(1-DAMP)*tbounce;
    !         // Reflect the position and velocity
    !         x[which] = 2*barrier-x[which];
    !         v[which] = -v[which];
    !         vh[which] = -vh[which];
    !         // Damp the velocities
    !         v[0] *= DAMP; vh[0] *= DAMP;
    !         v[1] *= DAMP; vh[1] *= DAMP;
    ! }

    !Coefficient of resitiution
    damp = 0.75

    !scale back the distance traveled based on time from collision
    tbounce = (sstate%x(i+which)-barrier) / sstate%v(i+which)

    sstate%x(i+0) = sstate%x(i+0) - sstate%v(i+0) * (1-damp)*tbounce
    sstate%x(i+1) = sstate%x(i+1) - sstate%v(i+1) * (1-damp)*tbounce

    ! reflect position and velocity
    sstate%x( i+which) = 2*barrier -sstate%x(i+which)
    sstate%v( i+which) =           -sstate%v(i+which)
    sstate%vh(i+which) =           -sstate%vh(i+which)


    ! damp velocities
    sstate%v(i+0)  = sstate%v(i+0)* damp
    sstate%vh(i+0) = sstate%vh(i+0)* damp
    sstate%v(i+1)  = sstate%v(i+1)* damp
    sstate%vh(i+1)  = sstate%vh(i+1)* damp


  end subroutine



    subroutine reflect_bc(sstate)
      use util
      type (systemstate) :: sstate

      integer :: i,n

      ! Boundaries of computational domain
      double precision :: xmin = 0.0d0
      double precision :: xmax = 1.0d0
      double precision :: ymin = 0.0d0
      double precision :: ymax = 1.0d0

      n = sstate%nParticles

      do i = 1,2*n,2 ! correct

        if (sstate%x(i) < xmin) then
          print *, "particle at left wall"
          call damp_reflect(i,0,xmin,sstate)
        end if

        if (sstate%x(i) > xmax) then
          print *, "particle at right wall"
          call damp_reflect(i,0,xmax,sstate)
        end if

        if (sstate%x(i+1) < ymin) then
          print *, "particle at bottom"
          call damp_reflect(i,1,ymin,sstate)
        end if

        if (sstate%x(i+1) > ymax) then
          print *, "particle at top"
          call damp_reflect(i,1,ymax,sstate)
        end if

      end do

    end subroutine




  subroutine compute_accel(sstate,params,ll,lc)
    use util
    use linkedlists

    type (systemstate)           :: sstate
    DOUBLE PRECISION, DIMENSION(9)  :: params
    integer, allocatable, dimension(:)        :: ll
    integer, allocatable, dimension(:,:)      :: lc
    integer                      :: n  !number of particles
    integer,dimension(2)         :: nmax
    integer                       :: i,j,no
    integer, dimension(4) :: ndx,ndy
    integer :: n1,n2,nx,ny
    double PRECISION :: h,rho0,k,mu,g,mass,h2,rcut
    double PRECISION :: c0,cp,cv
    integer          :: ncalcs  ! number of performed calculations (for performance testing)
    double PRECISION :: dx,dy,r2
    double precision :: rhoi,rhoj,q,u,w0,wp,wv,dvx,dvy

    n         = sstate%nParticles

    h         = params(3)
    rho0      = params(5)
    k         = params(6)
    mu        = params(7)
    g         = params(8)
    mass      = sstate%mass
    h2        = h*h

    rcut = params(9)             ! is 9th element in sim_param vector....
    nmax(1) = int(floor(1/rcut)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1/rcut)) ! maximum number of cells in each dimension

    !clearing ll and lc
    do i = 1,nmax(1)
      do j = 1,nmax(2)
        lc(i,j) = -1
      end do
    end do
    do i = 1,n
      ll(i) = -1
    end do


    ! start with gravity forces
    do i = 1,n
      sstate%a(2*i-1) = 0
      sstate%a(2*i-0) = (-1)*g
      ! print *, "sstate%a(2*i) = " ,sstate%a(2*i-0)
    end do

    ! constants for interaction term
    c0 = mass/pi/(h2*h2)
    cp = 15*k
    cv = -40*mu

    nCalcs = 0

    call setup_neighbour_list(sstate,params,ll,lc)
    call compute_density_with_ll(sstate,params,ll,lc)

    ndx = (/ 1,1,0,-1/)
    ndy = (/ 0,1,1,1 /)

    do i = 1,nmax(1)
      do j = 1,nmax(2)
        if (lc(i,j) /= -1) then
          n1 = lc(i,j)
          do while (n1 /= -1)
            n2 = ll(n1)
            rhoi = sstate%rho(n1)
            do while(n2 /= -1)
              nCalcs = nCalcs + 1
              dx = sstate%x(2*n2-1) - sstate%x(2*n1-1);
              dy = sstate%x(2*n2-0) - sstate%x(2*n1-0);
              r2 = dx*dx + dy*dy
              if(r2 < h2) then
                rhoj = sstate%rho(n2)
                q = sqrt(r2)/h
                u = 1-q
                w0 = c0*u/rhoi/rhoj
                wp = w0*cp*(rhoi+rhoj-2*rho0) * u/q
                wv = w0*cv
                dvx = sstate%v(2*n2-1) - sstate%v(2*n1-1)
                dvy = sstate%v(2*n2-0) - sstate%v(2*n1-0)

                sstate%a(2*n1-1) = sstate%a(2*n1-1) - (wp*dx + wv*dvx)
                sstate%a(2*n1-0) = sstate%a(2*n1-0) - (wp*dy + wv*dvy)
                sstate%a(2*n2-1) = sstate%a(2*n2-1) + (wp*dx + wv*dvx)
                sstate%a(2*n2-0) = sstate%a(2*n2-0) + (wp*dy + wv*dvy)
              end if
              n2 = ll(n2)
            end do

            !neighbor cells
            do no = 1,4
              nx = i + ndx(no)
              ny = j + ndy(no)

              ! boundaries
              if(nx < 0)         continue
              if(nx>nmax(1)-1)   continue
              if(ny<0)           continue
              if(ny>nmax(2)-1)   continue

              n2 = lc(nx,ny)

              do  while( n2 /= -1)
                nCalcs = nCalcs + 1
                dx = sstate%x(2*n2-1)-sstate%x(2*n1-1)
                dy = sstate%x(2*n2-0)-sstate%x(2*n1-0)
                r2 = dx*dx + dy*dy



                if (r2 < h2) then
                  rhoj = sstate%rho(n2)
                  q = sqrt(r2)/h
                  u = 1-q
                  !
                  ! print *, r2,h2,q,rhoi,rhoj,rho0

                  w0 = C0 * u/rhoi/rhoj
                  wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q
                  wv = w0 * Cv
                  dvx = sstate%v(2*n2-1) - sstate%v(2*n1-1)
                  dvy = sstate%v(2*n2-0) - sstate%v(2*n1-0)

                  ! print *, w0
                  ! print *, sstate%v

                  sstate%a(2*n1-1) = sstate%a(2*n1-1) - (wp*dx + wv*dvx)
                  sstate%a(2*n1-0) = sstate%a(2*n1-0) - (wp*dy + wv*dvy)
                  sstate%a(2*n2-1) = sstate%a(2*n1-1) + (wp*dx + wv*dvx)
                  sstate%a(2*n2-0) = sstate%a(2*n1-0) + (wp*dy + wv*dvy)


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





end module sphfunctions

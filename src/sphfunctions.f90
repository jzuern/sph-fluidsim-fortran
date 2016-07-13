module sphfunctions

  ! Implementation of functions and subroutines for the SPH method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 07/05/2016


implicit none

contains


  subroutine update_solid_particles_positions(sstate,params)

    ! updates solid particles position solely based on passed time

    use util
    type(systemstate)                           :: sstate !system state object
    type(sim_parameter)											    :: params
    double precision, dimension(3)              :: center = (/0.5d0, 0.4d0 , 0.5d0/) ! rotation center
    integer                                     :: particle
    double precision                            :: phi,phi_new, radius,dist_x,dist_y,dist_z

    do particle = 1 , sstate%nSolidParticles

      ! convert positions from cartesian to polar coordinates on order to make easy calculations of the new particle positions
      dist_x = sstate%x(3*particle-2) - center(1)
      dist_y = sstate%x(3*particle-1) - center(2)

      radius  =  SQRT (dist_x*dist_x + dist_y*dist_y )
      phi     =  ATAN2 (dist_x , dist_y)

      phi_new =  phi + params%dphi


      ! convert positions back to cartesian coordinates
      sstate%x(3*particle-2) = center(1) + radius * SIN(phi_new)
      sstate%x(3*particle-1) = center(2) + radius * COS(phi_new)


    end do

  end subroutine




  subroutine init_particles(sstate,params)

    ! initialize all particle positions

    use util
    type(systemstate)                           :: sstate !system state object
    type(sim_parameter)											    :: params

    call place_particles(sstate,params)
    call normalize_mass(sstate,params)

  end subroutine




  subroutine normalize_mass(sstate,params)

    ! In the initial configuration, we want the fluid to have the reference density.
    ! This is achieved by normalizing the total particle mass.

    use util

    type (systemstate)                          :: sstate
    type(sim_parameter)											    :: params
    double precision                            :: rho0
    double precision                            :: rho2s,rhos
    integer                                     :: idxstart, idxend

    rho0 = params%rho0


    ! get index range start and end corresponding to liquid particles
    idxstart = sstate%nSolidParticles + 1
    idxend   = sstate%nParticles

    call compute_density_initially(sstate,params) ! initial density calculation

    ! calculate sum of liquid particle densities
    rhos  = SUM(sstate%rho(idxstart : idxend))

    ! calculate sum of squared liquid particle densities
    rho2s = SUM(sstate%rho(idxstart : idxend) * sstate%rho(idxstart : idxend));


    sstate%mass = 1.d0
    sstate%mass = sstate%mass * (rho0*rhos / rho2s)

  end subroutine



  subroutine compute_density_with_ll(sstate, params, ll, lc)

    ! compute the density for each particle

    use util
    type (systemstate)                          :: sstate
    type(sim_parameter)											    :: params
    integer, allocatable, dimension(:)          :: ll
    integer, allocatable, dimension(:,:,:)      :: lc
    integer, dimension(26)                      :: ndx,ndy,ndz
    integer, dimension(3)                       :: nmax

    integer              :: nx,ny,nz
    integer              :: n1,n2,no
    integer              :: i,j,k
    double precision     :: mass,h,h2,h8,C,dx,dy,dz,r2,z,rho_ij

    h = params%h
    h2 = h*h
    h8 = h2*h2*h2*h2

    mass       = sstate%mass

    C = mass / pi / h8;

    ndx = (/1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1 /)
    ndy = (/1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1 /)
    ndz = (/1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1 /)

    nmax(1) = int(floor(1.d0/params%rcut_x)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/params%rcut_y)) ! maximum number of cells in y dimension
    nmax(3) = int(floor(1.d0/params%rcut_z)) ! maximum number of cells in y dimension


    !$omp parallel do private(n1,n2,dx,dy,dz,r2,no,nx,ny,nz,z,rho_ij)
    do i = 1,nmax(1) ! x-coordinate
      do j = 1,nmax(2) ! y-coordinate
      do k = 1,nmax(3) ! z-coordinate

        if (lc(i,j,k) /= -1) THEN
          n1 = lc(i,j,k)
          do while (n1 /= -1)
            n2 = ll(n1)
            sstate%rho(n1) = sstate%rho(n1) + 4*mass/Pi/h2

            do while(n2 /= -1)
              dx = sstate%x(3*n1-2) - sstate%x(3*n2-2)
              dy = sstate%x(3*n1-1) - sstate%x(3*n2-1)
              dz = sstate%x(3*n1-0) - sstate%x(3*n2-0)
              r2 = dx*dx + dy*dy + dz*dz
              z  = h2 - r2

              if (z > 0.0d0) then
                rho_ij = C*z*z*z
                sstate%rho(n1) = sstate%rho(n1) + rho_ij
                sstate%rho(n2) = sstate%rho(n2) + rho_ij
              end if

              n2 = ll(n2)
            end do

            ! Now the neighboring cells of cell i,j
            do no = 1,26
              nx = i+ndx(no)
              ny = j+ndy(no)
              nz = k+ndz(no)
              !boundary conditions
              if (nx <         1)  cycle
              if (nx > nmax(1)  )  cycle
              if (ny <         1)  cycle
              if (ny > nmax(2)  )  cycle
              if (nz <         1)  cycle
              if (nz > nmax(3)  )  cycle


              n2 = lc(nx,ny,nz)

              do while (n2 /= -1)
                dx = sstate%x(3*n1-2) - sstate%x(3*n2-2)
                dy = sstate%x(3*n1-1) - sstate%x(3*n2-1)
                dz = sstate%x(3*n1-0) - sstate%x(3*n2-0)
                r2 = dx*dx + dy*dy + dz*dz
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
    end do
    !  !$OMP END PARALLEL DO



  end subroutine

  subroutine compute_density_initially(sstate, params)


    use util
    type (systemstate)                       :: sstate
    type(sim_parameter)											 :: params
    integer, dimension(3)                    :: nmax
    integer                                  :: n,i,j
    double precision                         :: mass,h,h2,h8,C,dx,dy,dz,r2,z,rho_ij

    h = params%h
    h2 = h*h
    h8 = h2*h2*h2*h2

    n          = sstate%nParticles
    mass       = sstate%mass

    C = 4*mass / Pi / h8;


    nmax(1) = int(floor(1.d0/params%rcut_x)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/params%rcut_y)) ! maximum number of cells in y dimension
    nmax(3) = int(floor(1.d0/params%rcut_z)) ! maximum number of cells in y dimension

    do i = 1,n
      sstate%rho(i) = sstate%rho(i) + 4*mass/Pi/h2
      do j = i+1,n
        dx = sstate%x(3*i-2) - sstate%x(3*j-2)
        dy = sstate%x(3*i-1) - sstate%x(3*j-1)
        dz = sstate%x(3*i-0) - sstate%x(3*j-0)
        r2 = dx*dx + dy*dy + dz*dz
        z = h2-r2
        if (z > 0.d0) THEN
          rho_ij = C*z*z*z
          sstate%rho(i) = sstate%rho(i) + rho_ij
          sstate%rho(j) = sstate%rho(j) + rho_ij
        end if
      end do
    end do


  end subroutine




  subroutine place_particles(sstate,params)


    ! place the particles in the simulation domain.
    ! First, the number of particles based on their size and distance is determined.
    ! Then, we allocate the required fields for their position, velocity, ...
    ! After that we can actually save every initial particle position in those now allocated fields


    use util
    type (systemstate)                    :: sstate
    type(sim_parameter)										:: params
    double precision                      :: h,hh_liquid, hh_solid
    integer                               :: solidParticle,liquidParticle
    double precision                      :: x,y,z,rd
    integer                               :: nSolidParticles = 0, nLiquidParticles = 0



    h = params%h ! size of particles

    hh_liquid = 1.0d0*h ! distance of liquid particles from one another in the initial configuration
    hh_solid  = 1.0d0*h ! distance of solid particles from one another

    if (params%mill) THEN

      ! count solid particles
      x = 0.d0
      do while (x < 1.d0)
        y = 0.d0
        do while (y < 1.d0)
          z = 0.d0
          do while (z < 1.0d0)
            if (cross_indicator(x,y,z,params) /= 0) THEN
              nSolidParticles = nSolidParticles + 1
            end if
            z = z + hh_solid
          end do
          y = y + hh_solid
        end do
        x = x + hh_solid
      end do
    end if !mill


    sstate%nSolidParticles = nSolidParticles;

    ! count liquid particles
    x = 0.d0
    do while (x < 1.d0)
      y = 0.d0
      do while (y < 1.d0)
        z = 0.d0
        do while (z < 1.0d0)
          if (circ_indicator(x,y,z) /= 0 .AND. cross_indicator(x,y,z,params) == 0) THEN
            nLiquidParticles = nLiquidParticles + 1
          end if
          z = z + hh_liquid
        end do
        y = y + hh_liquid
      end do
      x = x + hh_liquid
    end do

    print *, "Number of liquid particles in simulation: " , nLiquidParticles
    print *, "Number of solid particles in simulation: " , nSolidParticles

    ! set number of liquid particles to counted value
    sstate%nLiquidParticles = nLiquidParticles

    ! set total number particles
    sstate%nParticles = sstate%nSolidParticles + nLiquidParticles

    call alloc_state(sstate,params)      ! allocate fields for x,v,vh,a

    solidParticle = 1 ! current particle

    if (params%mill) THEN

    ! place solid particles
      x = 0.d0
      do while (x < 1.d0)
        y = 0.d0
        do while (y < 1.d0)
          z = 0.0d0
          do while ( z < 1.0d0)
            if (cross_indicator(x,y,z,params) /= 0) THEN
              sstate%x(2*solidParticle-1) = x
              sstate%x(2*solidParticle-0) = y
              solidParticle = solidParticle + 1
            end if
            z = z + hh_solid
          end do
          y = y + hh_solid
        end do
        x = x + hh_solid
      end do
    end if

    liquidParticle = sstate%nSolidParticles + 1 ! current particle

    ! place liquid particles
    x = 0.d0
    do while (x < 1.d0)
      y = 0.d0
      do while (y < 1.d0)
        z = 0.d0
        do while (z < 1.d0)
          if (circ_indicator(x,y,z) /= 0 .AND. cross_indicator(x,y,z,params) == 0) THEN
            CALL RANDOM_NUMBER(rd)   ! random number between 0 and 1

            ! add some random noise to particle positions in order to prevent any
            ! symmetries to be preserved during simulation

            sstate%x(3*liquidParticle-2) = x + rd*0.0001d0
            sstate%x(3*liquidParticle-1) = y + rd*0.0001d0
            sstate%x(3*liquidParticle-0) = z + rd*0.0001d0

            sstate%v(3*liquidParticle-2) = 0.d0
            sstate%v(3*liquidParticle-1) = 0.d0
            sstate%v(3*liquidParticle-0) = 0.d0

            liquidParticle = liquidParticle + 1
          end if
          z = z + hh_liquid
        end do
        y = y + hh_liquid
      end do
      x = x + hh_liquid
    end do

  end subroutine


  subroutine damp_reflect(i,which, barrier, sstate)

    ! Auxiliary subroutine for the particle reflection that is called by reflect_bc()

    use util
    type (systemstate)  :: sstate
    integer             :: which,i
    double precision    :: barrier
    double precision    :: damp, tbounce

    ! ignore degenerate cases
    if (sstate%v(i-2+which) == 0.d0) then
      return
    end if

    !Coefficient of resitiution
    damp = 0.75d0

    !scale back the distance traveled based on time from collision
    tbounce = (sstate%x(i-2+which) - barrier) / sstate%v(i-2+which)

    sstate%x(i-2) = sstate%x(i-2) - sstate%v(i-2) * (1-damp)*tbounce
    sstate%x(i-1) = sstate%x(i-1) - sstate%v(i-1) * (1-damp)*tbounce
    sstate%x(i  ) = sstate%x(i  ) - sstate%v(i  ) * (1-damp)*tbounce

    ! reflect position and velocity
    sstate%x( i-2+which) = 2*barrier -sstate%x (i-2+which)
    sstate%v( i-2+which) =           -sstate%v( i-2+which)
    sstate%vh(i-2+which) =           -sstate%vh(i-2+which)

    ! damp velocities
    sstate%v (i-2)  = sstate%v (i-2)* damp
    sstate%vh(i-2)  = sstate%vh(i-2)* damp
    sstate%v (i-1)  = sstate%v (i-1)* damp
    sstate%vh(i-1)  = sstate%vh(i-1)* damp
    sstate%v( i  )  = sstate%v (i  )* damp
    sstate%vh(i  )  = sstate%vh(i  )* damp


  end subroutine



    subroutine reflect_bc(sstate)


      ! implementation of the reflection behavior on rigid boundaries


      use util
      type (systemstate) :: sstate

      integer :: i,first,last

      ! Boundaries of computational domain
      double precision :: xmin = 0.0d0
      double precision :: xmax = 1.0d0
      double precision :: ymin = 0.0d0
      double precision :: ymax = 1.0d0
      double precision :: zmin = 0.0d0
      double precision :: zmax = 1.0d0

      first = sstate%nSolidParticles+1
      last  = sstate%nParticles

      do i = first,last

        if (sstate%x(3*i - 2) < xmin) then
          ! print *, "particle at left wall"
          call damp_reflect(3*i,0,xmin,sstate)
        end if

        if (sstate%x(3*i - 2) > xmax) then
          ! print *, "particle at right wall"
          call damp_reflect(3*i,0,xmax,sstate)
        end if

        if (sstate%x(3*i - 1) < ymin) then
          ! print *, "particle at bottom"
          call damp_reflect(3*i,1,ymin,sstate)
        end if

        if (sstate%x(3*i - 1) > ymax) then
          ! print *, "particle at top"
          call damp_reflect(3*i,1,ymax,sstate)
        end if

        if (sstate%x(3*i - 0) < zmin) then
          call damp_reflect(3*i,2,zmin,sstate)
        end if

        if (sstate%x(3*i - 0) > zmax) then
          call damp_reflect(3*i,2,zmax,sstate)
        end if

      end do

    end subroutine




  subroutine compute_accel(sstate,params,ll,lc)

    ! compute the acceleration of the particles due to external forces and
    ! collisions with one another

    use util
    use linkedlists

    type (systemstate)                    :: sstate
    type(sim_parameter)										:: params
    integer, allocatable, dimension(:)    :: ll
    integer, allocatable, dimension(:,:,:)  :: lc
    integer                               :: n  !number of particles
    integer,dimension(3)                  :: nmax
    integer                               :: i,j,k,no
    integer,dimension(26)                  :: ndx,ndy,ndz
    integer                               :: n1,n2,nx,ny,nz
    double precision                      :: h,rho0,kk,mu,g,mass,h2
    double precision                      :: c0,cp,cv
    double precision                      :: dx,dy,dz,r2
    double precision                      :: rhoi,rhoj,q,u,w0,wp,wv,dvx,dvy,dvz
    double precision                      :: sp

    n         = sstate%nParticles
    h         = params%h
    rho0      = params%rho0
    kk        = params%k
    mu        = params%mu
    g         = params%g
    mass      = sstate%mass
    h2        = h*h

    nmax(1) = int(floor(1.d0/params%rcut_x)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/params%rcut_y)) ! maximum number of cells in y dimension
    nmax(3) = int(floor(1.d0/params%rcut_z)) ! maximum number of cells in y dimension

    !clearing ll and lc
    lc = -1
    ll = -1

    !set scale_parameter (experimental feature)
    sp = 1.0d0

    ! apply gravitational forces (=acceleration)
    sstate%a = 0.0d0
    sstate%a(2:3*n:3) = -g ! only apply to y-component of acceleration vector


    ! constants for interaction term
    c0 = mass/pi/(h2*h2)
    cp = 15.0d0*kk
    cv = -40.0d0*mu

    ! update neighbor list and density distribution
    call setup_neighbour_list(sstate,params,ll,lc)
    call compute_density_with_ll(sstate,params,ll,lc)

    ! print *, "in compute_accel 1:"
    ! print *, sstate%a

    ndx = (/1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1 /)
    ndy = (/1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0,-1,-1,-1,1,1,1,0,0,0,-1,-1,-1 /)
    ndz = (/1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1 /)

    ! Use openMP to parallelize cell access
    !$omp parallel do private(n1,n2,rhoi,dx,dy,dz,r2,rhoj,q,u,w0,wp,wv,dvx,dvy,dvz,no,nx,ny,nz)
    do i = 1,nmax(1)
      do j = 1,nmax(2)
      do k = 1,nmax(3)

        ! check for particles in cell (i,j)
        if (lc(i,j,k) /= -1) then
          n1 = lc(i,j,k)

          ! check for other particles in same cell as n1
          do while (n1 /= -1)
            n2 = ll(n1)
            rhoi = sstate%rho(n1)

            ! go through all particles in same cell as n1
            do while(n2 /= -1)
              dx = sstate%x(3*n2-2) - sstate%x(3*n1-2);
              dy = sstate%x(3*n2-1) - sstate%x(3*n1-1);
              dz = sstate%x(3*n2-0) - sstate%x(3*n1-0);
              r2 = dx*dx + dy*dy + dz*dz

              if(r2 < h2) then ! particles touching
                rhoj = sstate%rho(n2)
                q = sqrt(r2)/h
                u = 1-q
                w0 = c0*u/rhoi/rhoj
                wp = w0*cp*(rhoi+rhoj-2*rho0) * u/q
                wv = w0*cv
                dvx = sstate%v(3*n2-2) - sstate%v(3*n1-2)
                dvy = sstate%v(3*n2-1) - sstate%v(3*n1-1)
                dvz = sstate%v(3*n2-0) - sstate%v(3*n1-0)

                ! test if particle n1 is actually liquid particle.
                ! Then we can update acceleration accordingly
                if ( n1 > sstate%nSolidParticles) THEN
                  ! print *, sp,wp,dx,wv,dvx
                  sstate%a(3*n1-2) = sstate%a(3*n1-2) - sp*(wp*dx + wv*dvx)
                  sstate%a(3*n1-1) = sstate%a(3*n1-1) - sp*(wp*dy + wv*dvy)
                  sstate%a(3*n1-0) = sstate%a(3*n1-0) - sp*(wp*dz + wv*dvz)
                end if
                if (n2 > sstate%nSolidParticles) THEN
                  sstate%a(3*n2-2) = sstate%a(3*n2-2) + sp*(wp*dx + wv*dvx)
                  sstate%a(3*n2-1) = sstate%a(3*n2-1) + sp*(wp*dy + wv*dvy)
                  sstate%a(3*n1-0) = sstate%a(3*n1-0) + sp*(wp*dz + wv*dvz)
                end if


              end if
              n2 = ll(n2)
            end do

            !check for neighboring particles in neighbor cells as well
            do no = 1,26
              nx = i + ndx(no)
              ny = j + ndy(no)
              nz = k + ndz(no)
              !
              ! print *, nmax

              ! boundaries
              if(nx < 1)         cycle
              if(nx > nmax(1))   cycle
              if(ny < 1)         cycle
              if(ny > nmax(2))   cycle
              if(nz < 1)         cycle
              if(nz > nmax(3))   cycle

              ! print *, "nx ny nz:" ,nx,ny,nz

              n2 = lc(nx,ny,nz)

              do  while( n2 /= -1)

                dx = sstate%x(3*n2-2)-sstate%x(3*n1-2)
                dy = sstate%x(3*n2-1)-sstate%x(3*n1-1)
                dz = sstate%x(3*n2-0)-sstate%x(3*n1-0)
                r2 = dx*dx + dy*dy + dz*dz

                if (r2 < h2) then

                  rhoj = sstate%rho(n2)
                  q = sqrt(r2)/h
                  u = 1-q
                  w0 = c0 * u/rhoi/rhoj
                  wp = w0 * cp * (rhoi+rhoj-2*rho0) * u/q
                  wv = w0 * cv
                  dvx = sstate%v(3*n2-2) - sstate%v(3*n1-2)
                  dvy = sstate%v(3*n2-1) - sstate%v(3*n1-1)
                  dvz = sstate%v(3*n2-0) - sstate%v(3*n1-0)


                  ! test if particle n1 is actually liquid particle.
                  ! Then we can update acceleration accordingly
                  if ( n1 > sstate%nSolidParticles) THEN
                    sstate%a(3*n1-2) = sstate%a(3*n1-2) - sp*(wp*dx + wv*dvx)
                    sstate%a(3*n1-1) = sstate%a(3*n1-1) - sp*(wp*dy + wv*dvy)
                    sstate%a(3*n1-0) = sstate%a(3*n1-0) - sp*(wp*dz + wv*dvz)

                  end if
                  if (n2 > sstate%nSolidParticles) THEN
                    sstate%a(3*n2-2) = sstate%a(3*n2-2) + sp*(wp*dx + wv*dvx)
                    sstate%a(3*n2-1) = sstate%a(3*n2-1) + sp*(wp*dy + wv*dvy)
                    sstate%a(3*n2-0) = sstate%a(3*n2-0) + sp*(wp*dz + wv*dvz)

                  end if
                end if
                n2 = ll(n2)
              end do
            end do
            n1 = ll(n1)
          end do
        end if
      end do
      end do
    end do
    ! !$omp end parallel do

    ! print *, "in compute_accel 2:"
    ! print *, sstate%a
    !
    ! print *, " out compute_accel"

  end subroutine





end module sphfunctions

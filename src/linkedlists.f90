module linkedlists


  ! Implementation of functions and subroutines for linked lists bookkeeping

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 06/01/2016

  implicit none

contains


  subroutine setup_neighbour_list(sstate, params, ll,lc)

    use util
    type(systemstate)                                 :: sstate
    type(sim_parameter)											          :: params
    integer, allocatable, dimension(:)                :: ll
    integer, allocatable, dimension(:,:)              :: lc

    integer               :: i
    integer               :: ntot
    double precision      :: rcut
    integer, dimension(2) :: nmax, nidx

    ntot = sstate%nParticles   ! total number of particles
    rcut = params%rcut


    nmax(1) = int(floor(1.d0/rcut)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1.d0/rcut)) ! maximum number of cells in each dimension

    do i = 1,ntot
      nidx(1) = int(floor((sstate%x(2*i-1))/rcut)); !x coordinate
      nidx(1) = min(nidx(1),nmax(1))
      nidx(1) = max(nidx(1),1)

      nidx(2) = int(floor((sstate%x(2*i-0))/rcut)); !y coordinate
      nidx(2) = min(nidx(2),nmax(2))
      nidx(2) = max(nidx(2),1)

      ll(i) = lc(nidx(1),nidx(2))

      lc(nidx(1),nidx(2)) = i
    end do


  end subroutine


  subroutine print_neighour_list(sstate, params, ll,lc)

    ! print neighbor list for debugging purposes

    use util
    type(systemstate)              :: sstate
    type(sim_parameter)						 :: params
    integer, dimension(:)          :: ll
    integer, dimension(:,:)        :: lc

    integer               :: i,j,n
    integer               :: ntot
    double precision      :: rcut
    integer, dimension(2) :: nmax, nidx


    ntot = sstate%nParticles   ! total number of particles
    rcut = params%rcut


    nmax(1) = int(floor(1.d0/rcut)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/rcut)) ! maximum number of cells in y dimension

    do i = 1, nmax(1)
      do j = 1,nmax(2)
        if (lc(i,j) /= -1) THEN
          n = lc(i,j)
          print *, "cell ", i,j, ": "
          do while ( n /= -1)
             print*, "particle ", n , ",coordinates: " , sstate%x(2*n-1) , " " , sstate%x(2*n-0)
            n = ll(n)
          end do
          print *,
        else
          ! print *, " No particles in cell " , i , " " , j
        end if
      end do
    end do


  end subroutine


  subroutine init_ll(sstate, ll)

    use util
    type(systemstate)                     :: sstate
    integer, allocatable, dimension(:)    :: ll
    integer                               :: i
    integer                               :: ntot


    ntot  = sstate%nParticles
    allocate(ll(ntot))

    ll = -1 ! initialize ll with -1 (empty cell)

  end subroutine


  subroutine init_lc(sstate, params, lc)
    use util
    type(systemstate)                        :: sstate
    type(sim_parameter)									  	 :: params
    integer, allocatable, dimension(:,:)     :: lc
    double precision                         :: rcut
    integer, dimension(2)                    :: nmax
    integer                                  :: i,j
    integer                                  :: ntot

    ntot = sstate%nParticles
    rcut = params%rcut ! cutoff radius


    nmax(1) = int(floor(1.d0/rcut)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/rcut)) ! maximum number of cells in y dimension

    allocate(lc(nmax(1),nmax(2)))
    lc = -1 ! initialize lc with -1 (empty cell)


  end subroutine


end module linkedlists

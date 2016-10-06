module linkedlists


  ! Implementation of functions and subroutines for linked lists bookkeeping

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 07/05/2016

  implicit none

contains


  subroutine setup_neighbour_list(sstate, params, ll,lc)

    use util
    type(systemstate)                                 :: sstate
    type(sim_parameter)											          :: params
    integer, allocatable, dimension(:)                :: ll
    integer, allocatable, dimension(:,:,:)              :: lc

    integer               :: i
    integer               :: ntot
    double precision      :: rcut_x,rcut_y,rcut_z
    integer, dimension(3) :: nmax, nidx

    rcut_x = params%rcut_x
    rcut_y = params%rcut_y
    rcut_z = params%rcut_z


    nmax(1) = int(floor(1.d0/rcut_x)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1.d0/rcut_y)) ! maximum number of cells in each dimension
    nmax(3) = int(floor(1.d0/rcut_z)) ! maximum number of cells in each dimension

    do i = 1 , sstate%nParticles
      nidx(1) = int(floor((sstate%x(3*i-2))/rcut_x)); !x coordinate
      nidx(1) = min(nidx(1),nmax(1))
      nidx(1) = max(nidx(1),1)

      nidx(2) = int(floor((sstate%x(3*i-1))/rcut_y)); !y coordinate
      nidx(2) = min(nidx(2),nmax(2))
      nidx(2) = max(nidx(2),1)

      nidx(3) = int(floor((sstate%x(3*i-0))/rcut_z)); !z coordinate
      nidx(3) = min(nidx(3),nmax(3))
      nidx(3) = max(nidx(3),1)

      ll(i) = lc(nidx(1),nidx(2),nidx(3))

      lc(nidx(1),nidx(2),nidx(3)) = i
    end do


  end subroutine


  subroutine print_neighour_list(sstate, params, ll,lc)

    ! print neighbor list (for debugging purposes)

    use util
    type(systemstate)              :: sstate
    type(sim_parameter)						 :: params
    integer, dimension(:)          :: ll
    integer, dimension(:,:,:)        :: lc

    integer               :: i,j,k,n
    integer, dimension(3) :: nmax, nidx


    nmax(1) = int(floor(1.d0/params%rcut_x)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/params%rcut_y)) ! maximum number of cells in y dimension
    nmax(3) = int(floor(1.d0/params%rcut_z)) ! maximum number of cells in z dimension

    do i = 1, nmax(1)
      do j = 1,nmax(2)
      do k = 1,nmax(3)
        if (lc(i,j,k) /= -1) THEN
          n = lc(i,j,k)
          print *, "cell ", i,j,k, ": "
          do while ( n /= -1)
             print*, "particle ", n , ",coordinates: " , sstate%x(3*n-2) , " " , sstate%x(3*n-1), " " , sstate%x(3*n-0)
            n = ll(n)
          end do
          print *,
        else
          ! print *, " No particles in cell " , i , " " , j , " " , k
        end if
      end do
      end do
    end do

  end subroutine


  subroutine init_ll(sstate, ll)


    ! Initialize the linked list that keeps track of which particles sit in the same
    ! cell as the current particle

    use util
    type(systemstate)                     :: sstate
    integer, allocatable, dimension(:)    :: ll

    allocate(ll(sstate%nParticles))

    ll = -1 ! initialize ll to -1 (empty cell)

  end subroutine


  subroutine init_lc(sstate, params, lc)

    ! Initialize the field lc that contains the number of the first particle
    ! in a cell

    use util
    type(systemstate)                        :: sstate
    type(sim_parameter)									  	 :: params
    integer, allocatable, dimension(:,:,:)   :: lc
    integer, dimension(3)                    :: nmax

    nmax(1) = int(floor(1.d0/params%rcut_x)) ! maximum number of cells in x dimension
    nmax(2) = int(floor(1.d0/params%rcut_y)) ! maximum number of cells in y dimension
    nmax(3) = int(floor(1.d0/params%rcut_z)) ! maximum number of cells in y dimension

    allocate(lc(nmax(1),nmax(2),nmax(3)))
    lc = -1 ! initialize lc to -1 (empty cell)


  end subroutine


end module linkedlists

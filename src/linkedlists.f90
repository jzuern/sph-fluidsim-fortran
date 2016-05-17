module linkedlists


  ! Implementation of functions and subroutines for linked lists bookkeeping

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/16/2016



  implicit none

private
public :: setup_neighbour_list , print_neighour_list

contains


  subroutine setup_neighbour_list(sstate, params, ll,lc)
    use util
    type(systemstate)                                 :: sstate
    double precision,dimension(9),intent(in)          :: params
    integer, allocatable, dimension(:)                :: ll ! should already be allocated when passed (in theory...)
    integer, allocatable, dimension(:,:)              :: lc ! should already be allocated when passed

    integer               :: i
    integer               :: ntot
    double precision      :: rcut
    integer, dimension(2) :: nmax, nidx


    ntot = sstate%nParticles   ! total number of particles
    rcut = params(9)             ! is 9th element in sim_param vector....


    nmax(1) = int(floor(1/rcut)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1/rcut)) ! maximum number of cells in each dimension

    do i = 1,ntot
      nidx(1) = int((floor((sstate%x(2*i+0))/rcut))); !x coordinate
      nidx(1) = min(nidx(1),nmax(1)-1)
      nidx(1) = max(nidx(1),0)

      nidx(2) = int((floor((sstate%x(2*i+1))/rcut))); !y coordinate
      nidx(2) = min(nidx(2),nmax(2)-1)
      nidx(2) = max(nidx(2),0)

      ll(i) = lc(nidx(1),nidx(2))
      lc(nidx(1),nidx(2)) = i
    end do


  end subroutine


  subroutine print_neighour_list(sstate, params, ll,lc)
    ! print neighbor list for debugging purposes

    use util
    type(systemstate) :: sstate
    double precision, intent(in)   :: params
    integer, dimension(:)          :: ll
    integer, dimension(:,:)        :: lc




  end subroutine


end module linkedlists

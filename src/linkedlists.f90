module linkedlists


  ! Implementation of functions and subroutines for linked lists bookkeeping

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/19/2016

  ! To implement functions:
  ! - setup_neighbour_list
  ! - print_neighour_list
  ! - init_lc
  ! - init_ll



implicit none

private
public :: setup_neighbour_list , print_neighour_list , init_lc , init_ll

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
      nidx(1) = int((floor((sstate%x(2*i-1))/rcut))); !x coordinate
      nidx(1) = min(nidx(1),nmax(1)-1)
      nidx(1) = max(nidx(1),0)

      nidx(2) = int((floor((sstate%x(2*i-0))/rcut))); !y coordinate
      nidx(2) = min(nidx(2),nmax(2)-1)
      nidx(2) = max(nidx(2),0)

      ll(i) = lc(nidx(1),nidx(2))
      lc(nidx(1),nidx(2)) = i
    end do


  end subroutine


  subroutine print_neighour_list(sstate, params, ll,lc)
    ! print neighbor list for debugging purposes

    use util
    type(systemstate)              :: sstate
    double precision,dimension(9)  :: params
    integer, dimension(:)          :: ll
    integer, dimension(:,:)        :: lc

    integer               :: i,j,n
    integer               :: ntot
    double precision      :: rcut
    integer, dimension(2) :: nmax, nidx


    ntot = sstate%nParticles   ! total number of particles
    rcut = params(9)             ! is 9th element in sim_param vector....


    nmax(1) = int(floor(1/rcut)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1/rcut)) ! maximum number of cells in each dimension

    do i = 1, nmax(1)

      do j = 1,nmax(2)

        if (lc(i,j) /= -1) THEN
          n = lc(i,j)

          do while ( n /= -1)
            print *, n , ",coordinates: " , sstate%x(2*n-1) , " " , sstate%x(2*n-0)
            n = ll(n)
          end do
          print *, ! new line
        else
          print *, " No particles in cell " , i , " " , j
        end if
      end do
    end do

    ! for(int i=0; i<nmax[0]; i++) {
    !         for(int j=0; j<nmax[1]; j++) {
    !                 if(lc[i][j]!=-1) {
    !                         n=lc[i][j];
    !                         std::cout<<"cell i,j:"<<i<<","<<j<<std::endl;
    !                         while(n!=-1) {
    !                                 std::cout<< n << ", coordinates: " << state->x[2*n] << " " << state->x[2*n+1] << std::endl;
    !                                 n=ll[n];
    !                         }
    !                         std::cout << std::endl;
    !                 }
    !                 else{
    !                         std::cout<<"no particles in cell "<<i<<","<<j<<std::endl;
    !                 }
    !         }
    ! }



  end subroutine


  subroutine init_ll(sstate, params, ll) ! implement as function or as subroutine?

    use util
    type(systemstate)             :: sstate
    double precision, intent(in)   :: params
    integer, dimension(:)          :: ll
    integer                           :: i
    integer                        :: ntot


    ntot  = sstate%nParticles

    do i = 1,ntot
        ll(i) = -1
    end do


  end subroutine





  subroutine init_lc(sstate, params, lc)
    use util
    type(systemstate)                        :: sstate
    double precision                         :: params
    integer, dimension(:,:)                  :: lc
    double precision                         :: rcut
    integer, dimension(2)                    :: nmax
    integer                                  :: i,j
    integer                                  :: ntot

    ntot = sstate%nParticles
    rcut = params(9) ! cutoff radius

    nmax(1) = int(floor(1/rcut)) ! maximum number of cells in each dimension
    nmax(2) = int(floor(1/rcut)) ! maximum number of cells in each dimension

    do i = 1,nmax(1)
      do j = 1,nmax(2)
        lc(i,j) = -1      ! initialize lc with -1 (empty cell)
      end do
    end do


  end subroutine


end module linkedlists

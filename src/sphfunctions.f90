module sphfunctions

  ! Implementation of functions and subroutines for the SPH method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/16/2016



implicit none

contains
  subroutine reflect_bc(s)
    use util
    type (systemstate) :: s


  end subroutine



  subroutine compute_density_with_ll(s, params, ll, lc)
    use util
    type (systemstate)           :: s
    double precision, intent(in) :: params
    integer, dimension(:)        :: ll
    integer, dimension(:)        :: lc

  end subroutine


  subroutine compute_accel(s,params,ll,lc)
    use util
    type (systemstate)           :: s
    double precision, intent(in) :: params
    integer, dimension(:)        :: ll
    integer, dimension(:)        :: lc
  end subroutine












end module sphfunctions

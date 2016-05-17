module sphfunctions

  ! Implementation of functions and subroutines for the SPH method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/16/2016



implicit none


private

public :: reflect_bc, compute_density_with_ll, compute_accel,check_state

contains




    subroutine check_state(sstate)
      use util
      type(systemstate) :: sstate !system state object

      integer :: i
      double precision :: xi, yi

      do i = 1,sstate%nParticles
        xi = sstate%x(2*i + 0)
        yi = sstate%x(2*i + 1)
        ! assert statements unavailable in F90
        !         assert( xi >= 0 || xi <= 1 );
        !         assert( yi >= 0 || yi <= 1 );
      end do



    end subroutine




  subroutine reflect_bc(sstate)
    use util
    type (systemstate) :: sstate


  end subroutine



  subroutine compute_density_with_ll(s, params, ll, lc)
    use util
    type (systemstate)           :: s
    double precision, intent(in) :: params
    integer, dimension(:)        :: ll
    integer, dimension(:,:)        :: lc







  end subroutine


  subroutine compute_accel(s,params,ll,lc)
    use util
    type (systemstate)           :: s
    double precision, intent(in) :: params
    integer, dimension(:)        :: ll
    integer, dimension(:,:)        :: lc






  end subroutine












end module sphfunctions

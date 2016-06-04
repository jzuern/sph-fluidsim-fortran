module integrate

  ! Implementation of functions and subroutines for the numerical integration
  ! current method: leapfrog method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/24/2016

  implicit none


  contains

  subroutine leapfrog_start(state, dt)

    use util
    use sphfunctions

    type(systemstate)            :: state
    double precision,intent(in)  :: dt


    ! do i = 1,2*state%nParticles
    !   state%vh(i) = state%v(i) + state%a(i)*dt/2
    !   state%v(i)  = state%v(i) + state%a(i)*dt
    !   state%x(i)  = state%x(i) + state%vh(i)*dt
    ! end do

    state%vh = state%v  + (state%a   * dt/2)
    state%v  = state%v  + (state%a   * dt  )
    state%x  = state%x  + (state%vh  * dt  )

    call reflect_bc(state)

    return
  end


  subroutine leapfrog_step(state, dt)
    use util
    use sphfunctions   ! in order to call reflect_bc
    type (systemstate) :: state
    double precision   :: dt

    integer :: i
    !
    ! do i = 1,2*state%nParticles
    !   state%vh(i) = state%vh(i) + state%a(i)*dt
    !   state%v(i)  = state%vh(i) + state%a(i)*dt/2
    !   state%x(i)  = state%x(i)  + state%vh(i)*dt
    ! end do

    state%vh = state%vh + state%a*dt
    state%v  = state%vh + state%a*dt/2
    state%x  = state%x  + state%vh*dt


    call reflect_bc(state)

    return
  end





end module integrate

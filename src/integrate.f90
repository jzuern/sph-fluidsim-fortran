module integrate

  ! Implementation of functions and subroutines for the numerical integration
  ! current method: leapfrog method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 05/24/2016

  implicit none


  contains

  subroutine leapfrog_start(state, dt)
    use util
    use sphfunctions   ! in order to call reflect_bc
    type(systemstate) :: state
    double precision  :: dt

    print *, dt


    state%vh = state%vh + (state%a  * dt)
    state%v  = state%vh + (state%a  * dt/2)
    state%x  = state%x  + (state%vh * dt)


    ! const float*  a = s->a;
    ! float*  vh = s->vh;
    ! float*  v = s->v;
    ! float*  x = s->x;
    ! int n = s->n;
    !
    ! for (int i = 0; i < 2*n; ++i) vh[i] += a[i] * dt;
    ! for (int i = 0; i < 2*n; ++i) v[i] = vh[i] + a[i] * dt / 2;
    ! for (int i = 0; i < 2*n; ++i) x[i] += vh[i] * dt;
    !
    call reflect_bc(state)

    return
  end


  subroutine leapfrog_step(state, dt)
    use util
    use sphfunctions   ! in order to call reflect_bc
    type (systemstate) :: state
    double precision   :: dt

    state%vh = state%v + (state%a  * dt/2)
    state%v  = state%v + (state%a  * dt)
    state%x  = state%x + (state%vh * dt)

    ! const float* a = s->a;
    ! float*  vh = s->vh;
    ! float*  v = s->v;
    ! float*  x = s->x;
    ! int n = s->n;
    !
    ! for (int i = 0; i < 2*n; ++i) vh[i] = v[i] + a[i] * dt/2;
    ! for (int i = 0; i < 2*n; ++i) v[i] += a[i] * dt;
    ! for (int i = 0; i < 2*n; ++i) x[i] += vh[i] * dt;

    call reflect_bc(state)

    return
  end





end module integrate

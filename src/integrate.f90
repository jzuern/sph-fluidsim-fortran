module integrate

  ! Implementation of functions and subroutines for the numerical integration
  ! current method: leapfrog method

  !Created by Jannik Zuern on 05/16/2016
  !Last modified: 06/17/2016

  implicit none


  contains

  subroutine leapfrog_start(state,params)

    use util
    use sphfunctions

    type(systemstate)            :: state
    type(sim_parameter)					 :: params
    double precision             :: dt
    integer :: startidx, endidx

    ! only update liquid particle positions according to integration scheme
    startidx = 2*state%nSolidParticles+1
    endidx   = 2*state%nParticles

    dt = params%dt

    state%vh(startidx:endidx) = state%v(startidx:endidx)  + (state%a(startidx:endidx)   * dt/2)
    state%v(startidx:endidx)  = state%v(startidx:endidx)  + (state%a(startidx:endidx)   * dt  )
    state%x(startidx:endidx)  = state%x(startidx:endidx)  + (state%vh(startidx:endidx)  * dt  )

    call update_solid_particles_positions(state,params)
    call reflect_bc(state)

    return
  end


  subroutine leapfrog_step(state,params)
    use util
    use sphfunctions

    type (systemstate)             :: state
    type(sim_parameter)						 :: params
    double precision               :: dt
    integer :: startidx,endidx

    ! only update liquid particle positions according to integration scheme
    startidx = 2*state%nSolidParticles+1
    endidx   = 2*state%nParticles

    dt = params%dt

    state%vh(startidx:endidx) = state%vh(startidx:endidx) + state%a(startidx:endidx) * dt
    state%v(startidx:endidx)  = state%vh(startidx:endidx) + state%a(startidx:endidx) * dt/2
    state%x(startidx:endidx)  = state%x(startidx:endidx)  + state%vh(startidx:endidx)* dt

    call update_solid_particles_positions(state,params)
    call reflect_bc(state)

    return
  end


end module integrate

2D fluid simulation, based on Smoothed-Particle-Hydrodynamics (SPH).
Implementation: Linked Lists

Visualization:			 Gnuplot with Fortran interface ( http://www.math.yorku.ca/~akuznets/gnufor2/ )
				 Edited the source file gnufor2.f90 in order to improve the Gnuplot visualizations
Simulation parameters:    	 Command line argument to setup file setup.dat


// TODO:
- make parametertype from param_vector (otherwise very hard to read)
- separate src directory from .o / executable directory
- debug linked lists
- when to use interface?
- when to use public/private in module?
- when to use intent(in), intent(out), intent(inout) ?
- when to use integer, PARAMETER :: xyz
- ...

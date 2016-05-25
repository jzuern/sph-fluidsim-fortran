2D fluid simulation, based on Smoothed-Particle-Hydrodynamics (SPH).
Implementation: Linked Lists

Visualization:			 Gnuplot with Fortran interface ( http://www.math.yorku.ca/~akuznets/gnufor2/ )
				 Edited the source file gnufor2.f90 in order to modify the Gnuplot visualizations
Simulation parameters:    	 Command line argument to setup file setup.dat


// TODO:
SPECIFIC
- Debug:
	- cloud exploding when touching ground/when packed too much together. Spring forces too strong?!
		- possible explanations: variable not defined (= 0)

	- when program finishes: error in `./program': free(): invalid next size (normal): 0x00000000010f1e50





- make parametertype from param_vector (otherwise very hard to read)
- separate src directory from .o / executable directory


GENERAL
- when to use interface?
- when to use public/private in module?
- when to use intent(in), intent(out), intent(inout) ?
- when to use integer, PARAMETER :: xyz
- more usage of handy FORTRAN array wizardry?

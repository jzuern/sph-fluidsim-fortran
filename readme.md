Fortran 2D fluid simulation, based on Smoothed-Particle-Hydrodynamics (SPH).
Linked Lists

Visualization:			 
Gnuplot with Fortran interface ( http://www.math.yorku.ca/~akuznets/gnufor2/ )
				 Edited the source file gnufor2.f90 in order to modify the Gnuplot visualizations
Simulation parameters:    	 
Command line argument to setup file setup.dat


// BUILD

build the program with the provided MAKEFILE using

$ make


No external library dependencies

// USAGE

start the executable "program" with a simulation parameter file as the command line argument

$ ./program sim_parameter.dat

// EXPLANATION

This 2D SPH fluid simulation is loosely based on an introductary paper on this matter (http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf). I used the provided code fragments as an inspiration for my own implementation. Additionally I implemented a linked lists based approach in order to account for the quadratic time complexity when increasing the number of particles.


The code inside the source files is structured in a way as to reflect their role in the simulation.

-> gnufor2.f90: interface to call for Gnuplot instance inside FORTRAN program (credit to http://www.math.yorku.ca/~akuznets/gnufor2/)

-> integrate.f90: implementation of numerical integration technique

-> linkedlist.f90: implementation of linked list bookkeeping routines

-> sphfunctions:f90: implementation of particle interaction routines

-> util.f90: auxiliary functions with no other suitable place

-> program.f90: main program which starts simulation





// TODO:
SPECIFIC
- Debug:
	- cloud exploding when touching ground/when packed too much together. Spring forces too strong?!
		- possible explanations: variable not defined (= 0)

	- when program finishes: error in `./program': free(): invalid next size (normal): 0x00000000010f1e50




GENERAL
- when to use interface?
- when to use public/private in module?
- when to use intent(in), intent(out), intent(inout) ?
- when to use integer, PARAMETER :: xyz
- more usage of handy FORTRAN array wizardry?

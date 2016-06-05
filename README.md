# Fortran 2D SPH Fluid Simulation

Fortran 2D fluid simulation, based on Smoothed-Particle-Hydrodynamics (SPH) with linked lists implementation




## BUILD


build the program with the provided MAKEFILE using

```
$ make
```


## USAGE

start the executable "program" with a simulation parameter file as the command line argument
i.e.:
```
$ ./program sim_parameter.dat
```

## EXPLANATION

This 2D SPH fluid simulation is loosely based on an [introductory paper](http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf) on this matter . I used the provided code fragments as an inspiration for my own implementation. Additionally I implemented a linked lists based approach in order to account for the quadratic time complexity when increasing the number of particles.

The code inside the source files is structured in a way as to reflect their role in the simulation.

### gnufor2.f90: interface to call for Gnuplot instance inside FORTRAN program (credit to its )

### integrate.f90: implementation of numerical integration technique

### linkedlist.f90: implementation of linked list bookkeeping routines

### sphfunctions:f90: implementation of particle interaction routines

### util.f90: auxiliary functions with no other suitable place

### program.f90: main program which starts simulation


## Visualization:

Gnuplot with [Fortran interface](http://www.math.yorku.ca/~akuznets/gnufor2/) Edited the source file gnufor2.f90 in order to modify the Gnuplot visualizations.


## Simulation parameters:

Command line argument to setup file setup.dat

Exemplary good parameters are:
```
200        - Number of frames
100         - Number of steps per frame. Default: 100
0.004        - Size of particles (radius). Also: Distance of particles in initial configuration
0.0001      - Time step size Default: 1E-4
1000        - Reference density Default: 1000
100000        - Bulk modulus (Kompressionsmodul) Default: 1E3
0.01         - Viscosity Default: 0.1
9.81        - gravity strength Default: 9.81
0.01         - Cutoff radius (as fraction of total size of simulation grid for linked lists neighbor tracking
```
Just create a sim_setup.dat file in your working directory and copy those lines into it.

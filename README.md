# Fortran 2D SPH Fluid Simulation

Fortran 2D fluid simulation, based on Smoothed-Particle-Hydrodynamics (SPH) with an  implementation of a Linked Cells Algorithm.


## DEPENDENCIES

We need openMP to be installed on the system in order to link to the openMP libraries.
If you do not wish to parallelize, simply remove the `-lopenmp` linker flag.

## BUILD


build the program with the provided MAKEFILE using the make command:

```
$ make
```


## USAGE

start the executable `program` with a simulation parameter file as the command line argument
i.e.:
```
$ ./program sim_parameter.dat
```

## EXPLANATION

This 2D SPH fluid simulation is loosely based on an [introductory paper](http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf) on this matter . I used the provided code fragments as an inspiration for my own implementation. Additionally I implemented a linked lists based approach in order to account for the quadratic time complexity when increasing the number of particles.

The code within the source files is structured in a way as to reflect its role in the simulation.
Here is a list of the most important files and functions / subroutines:

- gnufor2.f90: interface to call for Gnuplot instance inside FORTRAN program
- integrate.f90: implementation of the leapfrog integration method
    ⋅⋅* leapfrog_start:
    ⋅⋅* leapfrog_start:
- linkedlist.f90: implementation of linked list bookkeeping routines
    ⋅⋅* print_neighour_list: prints neighbor list to std out. For debugging purposes only
    ⋅⋅* setup_neighbour_list: initializes and updates the field lc and vector ll that keep track of in which cells the particles sit
- sphfunctions:f90: implementation of particle interaction routines
    ⋅⋅*
- util.f90: auxiliary functions with no other suitable place
    ⋅⋅*
- program.f90: main program which starts simulation


## Visualization:

Gnuplot with [Fortran interface](http://www.math.yorku.ca/~akuznets/gnufor2/)

I slightly edited the source file gnufor2.f90 in order to modify the Gnuplot visualizations.


## Simulation parameters:

Command line argument to parameter file sim_parameter.dat
This file is read using FORTRAN's namelist routine

Exemplary parameters for water-like behavior are:
```
&SIMPARAMETER
nframes=200                 ! Number of frames
nSteps_per_frame=50         ! Number of steps per frame. Default: 100
h=0.01                      ! Size of particles (radius). Also: Distance of particles in initial configuration
dt=0.0001                   ! Time step size Default: 1E-4
rho0=1000.0                 ! Reference density Default: 1000
k=10000.0                   ! Bulk modulus Default: 1E5
mu=1.0                      ! Viscosity Default: 1.0
g=9.81                      ! gravity strength Default: 9.81
rcut_x=0.1                  ! Cutoff radius in x-direction for cells (as fraction of total size of simulation grid)
rcut_y=0.1                  ! Cutoff radius in y-direction for cells (as fraction of total size of simulation grid)
dphi=-0.001                 ! rotation speed of cross
/

```
Create a sim_parameter.dat file in your working directory and copy those lines into it.
(or use the provided file)


## TODO

- update and improve and elongate README.md
- add images to README.md ?

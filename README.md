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

This 2D SPH fluid simulation is loosely based on an [introductory paper](http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf) on this matter . I used the provided code fragments as an inspiration for my own implementation.

Additionally, I implemented a linked lists based approach in order to account for the quadratic time complexity when increasing the number of particles.



The simulation setup introduces a particle cloud (liquid blob) that is falling downwards due to gravitational forces. I also implemented a rotating structure (a watermill) in the center of the simulation domain. This structure causes the blob to disperse and break up.



The code within the source files is structured in a way as to reflect its role in the simulation.
Below is a list of source files:

- gnufor2.f90: interface to call for Gnuplot instance inside FORTRAN program
- integrate.f90: implementation of the leapfrog integration method
- linkedlist.f90: implementation of linked list bookkeeping routines
- sphfunctions:f90: implementation of particle interaction routines
- util.f90: auxiliary functions with no other suitable place
- program.f90: main program which starts simulation

Extensive code commentary should be enough to help understanding the tasks of each function and subroutine.


## Visualization:

Gnuplot with [Fortran interface](http://www.math.yorku.ca/~akuznets/gnufor2/)

I slightly edited the source file gnufor2.f90 in order to modify the Gnuplot visualizations.


## Simulation parameters:

Command line argument to parameter file sim_parameter.dat
This file is read using FORTRAN's namelist routine

Exemplary parameters for water-like behavior are:
```
&SIMPARAMETER
nframes=100                ! Number of frames
nSteps_per_frame=100       ! Number of steps per frame
h=0.010                    ! Size of particles (radius)
dt=0.0001                  ! Time step size
rho0=1000.0                ! Reference density
k=1000000.0                ! Bulk modulus
mu=0.1                     ! Viscosity
g=9.81                     ! gravity strength
rcut_x=0.1                 ! cell cutoff radius in x-direction (as fraction of total size of simulation grid)
rcut_y=0.1                 ! cell cutoff radius in y-direction (as fraction of total size of simulation grid)
mill=.FALSE.               ! Decide whether a watermill is in the computational domain or not
dphi=-0.0003               ! rotation speed of watermill (only applicable if mill == .true.)
/

```
Create a sim_parameter.dat file in your working directory and copy those lines into it.
(or use the provided file)

## Screenshots

![ScreenShot](https://raw.github.com/jzuern/sph-fluidsim-fortran/master/data/images/1.png)


![ScreenShot](https://raw.github.com/jzuern/sph-fluidsim-fortran/master/data/images/2.png)

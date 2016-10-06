# Fortran 3D SPH Fluid Simulation

Fortran 3D fluid simulation, based on Smoothed-Particle-Hydrodynamics (SPH) with an  implementation of a Linked Cells Algorithm.

## TODO:
- plot mill and particles in different colors
- add paraview plotting script

## DEPENDENCIES

We need openMP to be installed on the system in order to link to the openMP libraries.
If you do not wish to parallelize, simply remove the `-lopenmp` linker flag.

We use [gnuplotfortran](http://gnuplotfortran.sourceforge.net/) to invoke a [Gnuplot](http://www.gnuplot.info/) session from within Fortran.
Installation guidelines for these libraries can be found on the respective websites.

You can link to the gnuplotfortran libraries either by adding the `.so` files to your `LD_LIBRARY_PATH` or by copying the  `.so` and  `.mod` fles into the project folder.

The linker flags for gnuplotfortran are `-lfortranposix` and `-lgnuplotfortran`.

## BUILD


build the program with the provided MAKEFILE using the make command:

```
$ make
```


## USAGE

Start the simulation by executing `program` with a simulation parameter file as the command line argument
i.e.:
```
$ ./program sim_setup.dat
```

After the simulation has ended you may plot the particle data in all its glory by running the bash script:
```
$ ./plot_data3d.sh
```

## EXPLANATION

This 3D SPH fluid simulation is loosely based on an [introductory paper](http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf) on this matter. I used the provided code fragments as an inspiration for my own implementation.

Additionally, I implemented a linked lists based approach in order to account for the quadratic time complexity depending on the particle count.

The simulation setup introduces a particle cloud (liquid blob) that is falling downwards due to gravitational forces. I also implemented a rotating structure (call it a watermill) in the center of the simulation domain. This structure causes the blob to disperse and break up.


The code within the source files is structured in a way as to reflect its role in the simulation.
Below, we have the list of source files:

- gnufor2.f90: interface to call for a Gnuplot instance from within the FORTRAN program
- integrate.f90: implementation of the leapfrog integration method
- linkedlist.f90: implementation of linked list bookkeeping routines
- sphfunctions:f90: implementation of particle interaction routines
- util.f90: auxiliary functions
- program.f90: main program which starts simulation and contains the simulation loop

Additionally commentary is provided along the code within the source files.

In order to plot the particle data that has been written to file, I suggest using the provided bash script `plot_data3d.sh` which starts Gnuplot and plots the contents of the `data` folder. You might need to edit the script to make it suit your plotting needs. You also might want to make the bash script executable by typing

```
$ chmod +x plot_data3d.sh
```



## Visualization:

[Gnuplotfortran](http://gnuplotfortran.sourceforge.net/) for on-line visualization of the particle cloud for debugging and quick observation purposes.

If you wish to visualize to data that has been written to file during the simulation, you may use the provided bash script `plot_data3d.sh` to plot the particle data contained in those files.


## Simulation parameters:

Command line argument to parameter file sim_setup.dat
This file is read using FORTRAN's namelist routine

Exemplary parameters for water-like behavior are:
```
&SIMPARAMETER
nframes=100              ! Number of frames
nSteps_per_frame=50      ! Number of steps per frame. Default: 50
h=0.015                  ! Size of particles (radius). Also: Distance of particles in initial configuration
dt=1.0E-4                ! Time step size. Default: 1E-4
rho0=1.0E4               ! Reference density Default: 10000
k=1.0E6                  ! Bulk modulus. Default: 1000
mu=0.001                 ! Viscosity. Default: 0.1
g=9.81                   ! gravity strength. Default: 9.81
rcut_x=0.025             ! cell cutoff radius in x-direction (as fraction of total size of simulation grid)
rcut_y=0.025             ! cell cutoff radius in y-direction (as fraction of total size of simulation grid)
rcut_z=0.025             ! cell cutoff radius in z-direction (as fraction of total size of simulation grid)
mill=.TRUE.              ! Decide whether a watermill is in the computational domain or not
dphi=-0.0004             ! rotation speed of watermill (only applicable if mill == .true.)
/

```
Create a sim_setup.dat file in your working directory and copy those lines into it.
(or use the provided file)

## Screenshots

![ScreenShot](https://raw.github.com/jzuern/sph-fluidsim-fortran/master/data/images/1.jpg)
![ScreenShot](https://raw.github.com/jzuern/sph-fluidsim-fortran/master/data/images/2.jpg)
![ScreenShot](https://raw.github.com/jzuern/sph-fluidsim-fortran/master/data/images/3.jpg)

#ifndef SPH_FUNCS_H
#define SPH_FUNCS_H

#include <stdio.h>
#include <math.h>
#include <vector>
#include "assert.h"
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include <omp.h>
#include "gnuplot-iostream.h"

struct sim_param_t{
  const char* filename; // File name
   int nFrames;    // Number of frames
   int nStepsFrame;// Number of steps per frame
   float h;        // Size of particles
   float dt;       // Time step size
   float rho0;     // Reference density
   float k;        // Bulk modulus (?)
   float mu;       // Viscosity
   float g;       // gravity strength
   float rcut;    // cutoff radius for neighbors
};

struct sim_state_t {
  int n; // number of particles in simulation
  float mass; // particle mass (for each one the same)
  float*  rho; // densities
  float*  x; // positions
  float*  vh; // velocities (half-step) (from leapfrog integration)
  float*  v; // velocities
  float*  a; // accellerations
};

sim_param_t initialize_sim_param_t(std::vector<std::string> paramvector);

void reflect_bc(sim_state_t* s);

int box_indicator(float x, float y);
int circ_indicator(float x, float y);

void plotPoints(int frame,std::ofstream& out,int n, sim_state_t* s,sim_param_t params);

sim_state_t* alloc_state(int n); // alloc memory for sim_state_t
void free_state(sim_state_t* s); // free memory of sim_state:t

void compute_density_without_ll(sim_state_t* s, sim_param_t* params);
void compute_density_with_ll(sim_state_t* s, sim_param_t* params, int* ll, int **lc);
void compute_accel(sim_state_t* s, sim_param_t* params, int* ll, int **lc);

sim_state_t* init_particles(sim_param_t* param);
sim_state_t* place_particles(sim_param_t* param);

void normalize_mass(sim_state_t* s, sim_param_t* param);
void check_state(sim_state_t* s);
static void damp_reflect(int which, float barrier, float* x, float* v, float* vh);


#endif // SPH_FUNCS_H

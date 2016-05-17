#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "sph_funcs.h"

void leapfrog_step(sim_state_t* s, double dt);
void leapfrog_start(sim_state_t* s, double dt);


#endif

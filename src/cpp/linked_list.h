#ifndef _LINKED_LIST_H
#define _LINKED_LIST_H


#include <iostream>
#include <cmath>
#include "sph_funcs.h"


int** init_lc(sim_param_t param, sim_state_t* state);
int* init_ll(sim_param_t param, sim_state_t* state);
void setup_neighbour_list(sim_param_t* param, sim_state_t* state, int*ll, int**lc);
void print_neighour_list(sim_state_t* state,sim_param_t param,int *ll, int **lc);

#endif

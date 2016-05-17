#include "linked_list.h"


int* init_ll(sim_param_t param, sim_state_t* state){
        int ntot = state->n; // number of particles in simulation
        int* ll = new int[ntot];
        for(int i=0; i<ntot; i++) { // initialize with value -1 (empty)
                ll[i]=-1;
        }
        return ll;
}

int** init_lc(sim_param_t param, sim_state_t* state){
        std::cout << "in init_lc \n";
        int ntot = state->n; // number of particles in simulation
        float rcut = param.rcut;
        int nmax[2]; // number of rcut cells in each dimension (x and y)

        for(int i=0; i<2; i++) { // go through each dimension
                nmax[i]=int(floor(1/rcut));
        }

        int** lc = new int*[nmax[0]];
        lc[0]= new int[(nmax[0])*(nmax[1])];

        for(int i=1; i<nmax[0]; i++) {
                lc[i]=lc[0]+i*(nmax[1]);
        }

        for(int i=0; i<nmax[0]; i++) {
                for(int j=0; j<nmax[1]; j++) {
                        lc[i][j]=-1; // initialize with value -1 (empty)
                }
        }
        return lc;
}


void setup_neighbour_list(sim_param_t* param, sim_state_t* state, int *ll, int * *lc){
        // std::cout << "in setup_neighbour_list \n";

        int ntot = state->n; // number of particles in simulation
        float rcut = param->rcut; // rcut equal in each dimension
        int nmax[2];
        for(int i=0; i<2; i++) { // go through each dimension
                nmax[i]=int(floor(1/rcut));
        }

        int nidx[2];
        for(int i=0; i<ntot; ++i) { // map each particle to a cell
                nidx[0]=int(floor((state->x[2*i])/rcut)); // x coordinate
                nidx[0]=std::min(nidx[0],nmax[0]-1); // check for domain
                nidx[0]=std::max(nidx[0],0);

                nidx[1]=int(floor(state->x[2*i+1]/rcut)); // y coordinate
                nidx[1]=std::min(nidx[1],nmax[1]-1); // check for domain
                nidx[1]=std::max(nidx[1],0);

                ll[i] = lc[nidx[0]][nidx[1]];
                lc[nidx[0]][nidx[1]]=i;
        }
}


// for debugging
void print_neighour_list(sim_state_t* state,sim_param_t param,int *ll, int **lc){

        std::cout << "Printing neighbor list\n";
        int n;
        int ntot = state->n; // number of particles in simulation
        float rcut = param.rcut; // rcut equal in each dimension
        int nmax[2];
        for(int i=0; i<2; i++) { // go through each dimension
                nmax[i]=int(floor(1/rcut));
        }

        for(int i=0; i<nmax[0]; i++) {
                for(int j=0; j<nmax[1]; j++) {
                        if(lc[i][j]!=-1) {
                                n=lc[i][j];
                                std::cout<<"cell i,j:"<<i<<","<<j<<std::endl;
                                while(n!=-1) {
                                        std::cout<< n << ", coordinates: " << state->x[2*n] << " " << state->x[2*n+1] << std::endl;
                                        n=ll[n];
                                }
                                std::cout << std::endl;
                        }
                        else{
                                std::cout<<"no particles in cell "<<i<<","<<j<<std::endl;
                        }
                }
        }
}

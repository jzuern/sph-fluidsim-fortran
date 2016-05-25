#include "sph_funcs.h"
#include "linked_list.h"

// #include "gnuplot.h" // must be included here..
// Gnuplot gpl;

#include <boost/tuple/tuple.hpp>

// This must be defined before the first time that "gnuplot-iostream.h" is included.
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"
Gnuplot gp;

sim_param_t initialize_sim_param_t(std::vector<std::string> paramvector){
        sim_param_t params;
        params.nFrames =        std::stoi(paramvector[0]);// Number of frames
        params.nStepsFrame =    std::stoi(paramvector[1]);// Number of steps per frame. Default: 100
        params.h =              std::stof(paramvector[2]);  // Size of particles Default: 5E-03
        params.dt =             std::stof(paramvector[3]); // Time step size Default: 1E-4
        params.rho0 =           std::stof(paramvector[4]);// Reference density Default: 1E3
        params.k  =             std::stof(paramvector[5]);  // Bulk modulus (Kompressionsmodul) Default: 1E3
        params.mu =             std::stof(paramvector[6]); // Viscosity Default: 0.1
        params.g =              std::stof(paramvector[7]); // gravity strength Default: 9.81
        params.rcut =           std::stof(paramvector[8]);// cutoff radius. Default: 0.01
        return params;
}



void plotPoints(int frame,std::ofstream& out, int n, sim_state_t* s,sim_param_t params){



        float* x = s->x;
        std::vector<boost::tuple<float,float> > points;
        for(int i = 0; i < n; i++) {
                points.push_back(boost::make_tuple(x[2*i+0], x[2*i+1]));
        }

        gp << "set xrange [0:1]\n";
        gp << "set yrange [0:1]\n";
        gp << "plot '-' with points pointtype 7\n";
        gp.send1d(points);

}


void compute_density_without_ll(sim_state_t* s, sim_param_t* params){
        int n = s->n;
        float*   rho = s->rho;
        float*   x = s->x;
        float mass = s->mass;
        float h = params->h;
        float h2 = h*h;
        float h8 = (h2*h2)*(h2*h2);
        float C = 4*s->mass / M_PI / h8;
        memset(rho, 0, n*sizeof(float)); // sets values of rho to (double)0

        for (int i = 0; i < n; i++) {
                rho[i] += 4*mass/M_PI/h2;
                for (int j = i+1; j < n; j++) {
                        float dx = x[2*i + 0] - x[2*j + 0];
                        float dy = x[2*i + 1] - x[2*j + 1];
                        float r2 = dx*dx + dy*dy;
                        float z = h2 - r2;
                        if (z > 0) {
                                float rho_ij = C*z*z*z;
                                rho[i] += rho_ij;
                                rho[j] += rho_ij;
                        }
                }
        }
}

void compute_density_with_ll(sim_state_t* s, sim_param_t* params, int* ll, int **lc){
        int n = s->n;
        float*   rho = s->rho;
        float*   x = s->x;
        float mass = s->mass;
        float h = params->h;
        float h2 = h*h;
        float h8 = (h2*h2)*(h2*h2);
        float C = 4*s->mass / M_PI / h8;
        memset(rho, 0, n*sizeof(float)); // sets values of rho to (double)0

        int nmax[2];
        int n1,n2;
        float rcut = params->rcut; // rcut equal in each dimension
        for(int i=0; i<2; i++) { // go through each dimension
                nmax[i]=int(floor(1/rcut));
        }


        int ndx[]={1,1,0,-1};
        int ndy[]={0,1,1,1};
        int nx,ny;

        for(int i=0; i<nmax[0]; i++) {
                for(int j=0; j<nmax[1]; j++) {
                        if(lc[i][j]!=-1) {
                                n1=lc[i][j];
                                while(n1 != -1) {
                                        n2=ll[n1];
                                        rho[n1] += 4*mass/M_PI/h2;
                                        while(n2!=-1) {
                                                float dx = x[2*n1+0]-x[2*n2+0];
                                                float dy = x[2*n1+1]-x[2*n2+1];
                                                float r2 = dx*dx + dy*dy;
                                                float z = h2 - r2;
                                                if (z > 0) {

                                                        float rho_ij = C*z*z*z;

                                                        // std::cout << "x[2*n1+0] = " << x[2*n1+0] << " x[2*n1+1] = " << x[2*n1+1] << "x[2*n2+0] = " << x[2*n2+0] << "x[2*n2+1] = " << x[2*n2+1] << std::endl;
                                                        rho[n1] += rho_ij;
                                                        rho[n2] += rho_ij;
                                                }
                                                // std::cout << "rho[i] = " << rho[n1] << ", rho[j] = " << rho[n2] << std::endl;
                                                n2=ll[n2];
                                        }
                                        // Neighbpr cells
                                        for(int no=0; no<4; no++) {
                                                nx=i+ndx[no];
                                                ny=j+ndy[no];

                                                // Randbedigungen:
                                                if(nx<0)           continue; // to end of loop...
                                                if(nx>nmax[0]-1)   continue;
                                                if(ny<0)           continue;
                                                if(ny>nmax[1]-1)   continue;
                                                // std::cout << "    checking neighbor cell " << nx << " " << ny << std::endl;
                                                n2=lc[nx][ny];

                                                while(n2!=-1) {
                                                        float dx = x[2*n1+0]-x[2*n2+0];
                                                        float dy = x[2*n1+1]-x[2*n2+1];
                                                        float r2 = dx*dx + dy*dy;
                                                        float z = h2 - r2;
                                                        if (z > 0) {
                                                                float rho_ij = C*z*z*z;
                                                                rho[n1] += rho_ij;
                                                                rho[n2] += rho_ij;
                                                        }
                                                        n2=ll[n2];
                                                }
                                        }
                                        n1=ll[n1];
                                }
                        }
                }
        }
}


void compute_accel(sim_state_t* state, sim_param_t* params, int* ll, int **lc){

        // Unpack basic parameters
        const float h= params->h;
        const float rho0 = params->rho0;
        const float k= params->k;
        const float mu= params->mu;
        const float g= params->g;
        const float mass = state->mass;
        const float h2= h*h;

        // Unpack system state
        const float*   rho = state->rho;
        const float*   x = state->x;
        const float*   v = state->v;
        float*         a = state->a;
        int n = state->n;


        int nmax[2];

                float rcut = params->rcut; // rcut equal in each dimension
                for(int i=0; i<2; i++) { // go through each dimension
                        nmax[i]=int(floor(1/rcut));
                }

                // clearing ll and lc:
                for(int i=0; i<nmax[0]; i++) {
                        for(int j=0; j<nmax[1]; j++) {
                                lc[i][j]=-1; // initialize with value -1 (empty)
                        }
                }
                for(int i=0; i<n; i++) { // initialize with value -1 (empty)
                        ll[i]=-1;
                }


        // Start with gravity and surface forces
        for (int i = 0; i < n; ++i) {
                a[2*i+0] = 0;
                a[2*i+1] = -g;
        }

        // Constants for interaction term
        float C0 = mass / M_PI / ( (h2)*(h2) );
        float Cp = 15*k;
        float Cv = -40*mu;

        int nCalcs = 0;


                // setting up updated neighbor list
                setup_neighbour_list(params,state,ll,lc); // update linked lists for new particle positions

                compute_density_with_ll(state,params,ll,lc);

                int ndx[]={1,1,0,-1};
                int ndy[]={0,1,1,1};
                int n1,n2,nx,ny;

                for(int i=0; i<nmax[0]; i++) {
                        // std::cout << "i = " << i << "\n";
                        for(int j=0; j<nmax[1]; j++) {
                                // std::cout << "j = " << j << "\n";
                                if(lc[i][j]!=-1) {
                                        n1=lc[i][j];
                                        while(n1 != -1) {
                                                n2=ll[n1];
                                                const float rhoi = rho[n1];
                                                while(n2!=-1) {
                                                        nCalcs += 1;
                                                        float dx = x[2*n2+0]-x[2*n1+0];
                                                        float dy = x[2*n2+1]-x[2*n1+1];
                                                        float r2 = dx*dx + dy*dy;
                                                        if (r2 < h2) {
                                                                const float rhoj = rho[n2];
                                                                float q = sqrt(r2)/h;
                                                                float u = 1-q;
                                                                float w0 = C0 * u/rhoi/rhoj;
                                                                float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                                                                float wv = w0 * Cv;
                                                                float dvx = v[2*n2+0]-v[2*n1+0];
                                                                float dvy = v[2*n2+1]-v[2*n1+1];
                                                                a[2*n1+0] -= (wp*dx + wv*dvx);
                                                                a[2*n1+1] -= (wp*dy + wv*dvy);
                                                                a[2*n2+0] += (wp*dx + wv*dvx);// muss drin bleiben
                                                                a[2*n2+1] += (wp*dy + wv*dvy);// muss drin bleiben
                                                        }
                                                        n2=ll[n2];
                                                }
                                                // Neighbor cells
                                                for(int no=0; no<4; no++) {
                                                        nx=i+ndx[no];
                                                        ny=j+ndy[no];

                                                        // boundaries:
                                                        if(nx<0)           continue; // to end of loop...
                                                        if(nx>nmax[0]-1)   continue;
                                                        if(ny<0)           continue;
                                                        if(ny>nmax[1]-1)   continue;
                                                        n2=lc[nx][ny];

                                                        while(n2!=-1) {
                                                                nCalcs += 1;
                                                                float dx = x[2*n2+0]-x[2*n1+0];
                                                                float dy = x[2*n2+1]-x[2*n1+1];
                                                                float r2 = dx*dx + dy*dy;
                                                                if (r2 < h2) {
                                                                        const float rhoj = rho[n2];
                                                                        float q = sqrt(r2)/h;
                                                                        float u = 1-q;
                                                                        float w0 = C0 * u/rhoi/rhoj;
                                                                        float wp = w0 * Cp * (rhoi+rhoj-2*rho0) * u/q;
                                                                        float wv = w0 * Cv;
                                                                        float dvx = v[2*n2+0]-v[2*n1+0];
                                                                        float dvy = v[2*n2+1]-v[2*n1+1];
                                                                        a[2*n1+0] -= (wp*dx + wv*dvx);
                                                                        a[2*n1+1] -= (wp*dy + wv*dvy);
                                                                        a[2*n2+0] += (wp*dx + wv*dvx);
                                                                        a[2*n2+1] += (wp*dy + wv*dvy);
                                                                }
                                                                n2=ll[n2];
                                                        }
                                                }
                                                n1=ll[n1];
                                        }
                                }
                        }
                }
        }
}


void free_state(sim_state_t* s){
        delete[] s->x;
        delete[] s->rho;
        delete[] s->a;
        delete[] s->v;
        delete[] s->vh;
        delete s;
}

int box_indicator(float x, float y){
        return (x > 0.3) && (y > 0.3) && (x < 0.7) && (y < 0.7);
}

int circ_indicator(float x, float y){
        float dx = (x-0.5);
        float dy = (y-0.6);
        float r2 = dx*dx + dy*dy;
        return (r2 > 0.1*0.1 && r2 < 0.3*0.3);
}


void normalize_mass(sim_state_t* s, sim_param_t* param){
        s->mass = 1;
        compute_density_without_ll(s, param);

        float rho0 = param->rho0;
        float rho2s = 0;
        float rhos = 0;
        for (int i = 0; i < s->n; ++i) {
                rho2s += (s->rho[i])*(s->rho[i]);
                rhos += s->rho[i];
        }
        s->mass *= ( rho0*rhos / rho2s );
}


sim_state_t* init_particles(sim_param_t* param){
        sim_state_t* s = place_particles(param);
        normalize_mass(s, param);
        return s;
}

sim_state_t* place_particles(sim_param_t* param){
        float h = param->h;
        float hh = h/1.0;// Initialize points closer than radius to each other so that
                         //realistic behaviour is obtained

        // Count mesh points that fall in indicated region.
        int count = 0;
        for (float x = 0; x < 1; x += hh) {
                for (float y = 0; y < 1; y += hh) {
                        count += circ_indicator(x,y);
                }
        }

        // Populate the particle data structure
        sim_state_t* s = alloc_state(count);

        // debug segmentation fault
        std::cout << " Number of particles in Simulation: " << s->n << std::endl;

        float rd;
        int p = 0;
        for (float x = 0; x < 1; x += hh) {
                for (float y = 0; y < 1; y += hh) {
                        if (circ_indicator(x,y)) {
                                rd = float(std::rand())/100000000000;
                                s->x[2*p+0] = x;
                                s->x[2*p+1] = y;
                                s->v[2*p+0] = rd;
                                s->v[2*p+1] = rd;
                                ++p;
                        }
                }
        }
        return s;
}

sim_state_t* alloc_state(int n){ // alloc memory for sim_state_t

        sim_state_t* state = new sim_state_t;
        float* rho = new float[n];
        float* v = new float[2*n];
        float* vh = new float[2*n];
        float* x = new float[2*n]; // 2 Dimensions (x,z)
        float* a = new float[2*n];

        state->n = n;
        state->rho = rho;
        state->v = v;
        state->vh = vh;
        state->x = x;
        state->a = a;
        state->mass = 0;

        return state;
}

void check_state(sim_state_t* s){ //check whether particles are contained within domain (just for debugging)
        for (int i = 0; i < s->n; ++i) {
                float xi = s->x[2*i+0];
                float yi = s->x[2*i+1];
                assert( xi >= 0 || xi <= 1 );
                assert( yi >= 0 || yi <= 1 );
        }
}


static void damp_reflect(int which, float barrier, float* x, float* v, float* vh) {
        // Coefficient of resitiution
        const float DAMP = 0.75;
        // Ignore degenerate cases
        if (v[which] == 0) return;
        // Scale back the distance traveled based on time from collision
        float tbounce = (x[which]-barrier)/v[which];
        x[0] -= v[0]*(1-DAMP)*tbounce;
        x[1] -= v[1]*(1-DAMP)*tbounce;
        // Reflect the position and velocity
        x[which] = 2*barrier-x[which];
        v[which] = -v[which];
        vh[which] = -vh[which];
        // Damp the velocities
        v[0] *= DAMP; vh[0] *= DAMP;
        v[1] *= DAMP; vh[1] *= DAMP;
}


void reflect_bc(sim_state_t* s){
        // Boundaries of the computational domain
        const float XMIN = 0.0;
        const float XMAX = 1.0;
        const float YMIN = 0.0;
        const float YMAX = 1.0;

        float*  vh = s->vh;
        float*  v = s->v;
        float*  x = s->x;
        int n = s->n;
        for (int i = 0; i < n; ++i, x += 2, v += 2, vh += 2) {
                if (x[0] < XMIN) damp_reflect(0, XMIN, x, v, vh);
                if (x[0] > XMAX) damp_reflect(0, XMAX, x, v, vh);
                if (x[1] < YMIN) damp_reflect(1, YMIN, x, v, vh);
                if (x[1] > YMAX) damp_reflect(1, YMAX, x, v, vh);
        }
}

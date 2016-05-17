#include "integrate.h"

void leapfrog_start(sim_state_t* s, double dt){
        const float* a = s->a;
        float*  vh = s->vh;
        float*  v = s->v;
        float*  x = s->x;
        int n = s->n;

        for (int i = 0; i < 2*n; ++i) vh[i] = v[i] + a[i] * dt/2;
        for (int i = 0; i < 2*n; ++i) v[i] += a[i] * dt;
        for (int i = 0; i < 2*n; ++i) x[i] += vh[i] * dt;

        reflect_bc(s);
}

void leapfrog_step(sim_state_t* s, double dt){
        const float*  a = s->a;
        float*  vh = s->vh;
        float*  v = s->v;
        float*  x = s->x;
        int n = s->n;

        for (int i = 0; i < 2*n; ++i) vh[i] += a[i] * dt;
        for (int i = 0; i < 2*n; ++i) v[i] = vh[i] + a[i] * dt / 2;
        for (int i = 0; i < 2*n; ++i) x[i] += vh[i] * dt;

        reflect_bc(s);
}

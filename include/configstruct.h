#ifndef CONFIGSTRUCT_H
#define CONFIGSTRUCT_H

#include "include/solvers.h"

typedef struct ConfigStruct
{
    double y0[6];
    double y_target[6];
    void (*propulsion_model)(double t, double y[6], double angles[2], double acceleration[3]);
    AdaptiveSolver *solver;
    void (*steering_law)(double t, double y[6], ConfigStruct *config, double angles[2]);
    double t_span[2];
    double ode_rel_tol;
    double guidance_tol;
    double guidance_weights[5];
    double penalty_param;
    double min_pe;
    double penalty_weight;
} ConfigStruct;

#endif // CONFIGSTRUCT_H
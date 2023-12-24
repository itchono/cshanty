#ifndef TYPES_H
#define TYPES_H

typedef struct RKSolution
{
    double *t;       // timesteps
    double (*y)[6];  // function value
    int n;           // number of steps taken
    int n_fev;       // number of function evaluations
    int n_step_fail; // number of failed steps
} RKSolution;

typedef struct ConfigStruct ConfigStruct;

typedef RKSolution(*(*ODESolver)(void(double, double *, double *, bool *, ConfigStruct *), double, double, ConfigStruct *));
typedef void (*SteeringLaw)(double, double *, ConfigStruct *, double[2]);
typedef void (*PropulsionModel)(double, double *, double[2], double[3]);

struct ConfigStruct
{
    double y0[6];
    double y_target[6];
    PropulsionModel propulsion_model;
    ODESolver solver;
    SteeringLaw steering_law;
    double t_span[2];
    double ode_rel_tol;
    double ode_h0;
    double guidance_tol;
    double guidance_weights[5];
    double penalty_param;
    double min_pe;
    double penalty_weight;
};

#endif
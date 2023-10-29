// Numerical Runge Kutta Solvers

#ifndef SOLVERS_H
#define SOLVERS_H

#define VEC_SIZE 6 // dimension of vector (kept static at 6 for orbit prop)

// Big macro for defining general RK method
#define RK_METHOD(weights_a, weights_b, weights_bh,                             \
                  weights_c, SOLN_ORD, NUM_STAGES)                              \
    if (tf <= t0)                                                               \
        return NULL; /* error case */                                           \
                                                                                \
    /* Initial estimate for n = number of rows in solution */                   \
    int n = (int)((tf - t0) / h0);                                              \
    int step = 0;                                                               \
                                                                                \
    /* Allocate result struct */                                                \
    RKSolution *result = (RKSolution *)malloc(sizeof(RKSolution));              \
    result->y = (double(*)[6])malloc(n * sizeof(double[VEC_SIZE]));             \
    result->t = (double(*))malloc(n * sizeof(double));                          \
                                                                                \
    /* Prepare Step variables */                                                \
    double h = h0;                  /* initial stepsize */                      \
    double t = t0;                  /* intial time */                           \
    double k[NUM_STAGES][VEC_SIZE]; /* k values; derivatives at each stage */   \
    double y_stage[VEC_SIZE];       /* used to evaluate derivs at each stage */ \
    double y_curr[VEC_SIZE];        /* solution at beginning of current step */ \
                                                                                \
    /* Init y0, t0 (step 0) */                                                  \
    for (int j = 0; j < VEC_SIZE; j++)                                          \
    {                                                                           \
        result->y[0][j] = y0[j];                                                \
        y_curr[j] = y0[j];                                                      \
    }                                                                           \
    result->t[0] = t;                                                           \
                                                                                \
    /* Create Step Variables */                                                 \
                                                                                \
    double err; /* magnitude of error, used for error estimate */               \
                                                                                \
    /* zero-out k */                                                            \
    for (int stage = 0; stage < NUM_STAGES; stage++)                            \
        for (int j = 0; j < VEC_SIZE; j++)                                      \
            k[stage][j] = 0;                                                    \
                                                                                \
    /* INTEGRATION STEPS */                                                     \
    while (t - h <= tf)                                                         \
    {                                                                           \
        /* Evaluate derivative at stage 0 */                                    \
        f(t, y_curr, k[0]);                                                     \
                                                                                \
        /*Evaluate derivatives at stages 1-15 */                                \
        for (int stage = 1; stage < NUM_STAGES; stage++)                        \
        {                                                                       \
            /* Prepare input vector */                                          \
            for (int j = 0; j < VEC_SIZE; j++)                                  \
            {                                                                   \
                y_stage[j] = y_curr[j]; /* take current sol */                  \
                for (int w = 0; w < stage; w++)                                 \
                    /* Add previous steps */                                    \
                    y_stage[j] += h * k[w][j] * weights_a[stage][w];            \
            }                                                                   \
            /* evaluate next k */                                               \
            f(t + h * weights_c[stage], y_stage, k[stage]);                     \
        }                                                                       \
                                                                                \
        /* Calculate solution estimates */                                      \
        err = 0;                                                                \
        double sol_lower_ord[VEC_SIZE];                                         \
        double sol_higher_ord[VEC_SIZE];                                        \
                                                                                \
        for (int j = 0; j < VEC_SIZE; j++)                                      \
        {                                                                       \
            sol_higher_ord[j] = y_curr[j];                                      \
            sol_lower_ord[j] = y_curr[j];                                       \
            for (int stage = 0; stage < NUM_STAGES; stage++)                    \
            {                                                                   \
                /* use b and b_h weights to calculate error */                  \
                sol_higher_ord[j] += h * weights_b[stage] * k[stage][j];        \
                sol_lower_ord[j] += h * weights_bh[stage] * k[stage][j];        \
            }                                                                   \
            err += pow(sol_higher_ord[j] - sol_lower_ord[j], 2);                \
        }                                                                       \
        err = sqrt(err);                                                        \
                                                                                \
        if (err < tol)                                                          \
        {                                                                       \
            /* Accept step */                                                   \
                                                                                \
            /* Check array size and increase size if needed */                  \
            if ((step++) >= n)                                                  \
            {                                                                   \
                n *= 2; /* double size of array and reallocate memory */        \
                result->y = (double(*)[VEC_SIZE])realloc(                       \
                    result->y, n * sizeof(double[VEC_SIZE]));                   \
                result->t = (double(*))realloc(result->t, n * sizeof(double));  \
            }                                                                   \
                                                                                \
            /* Update t, y_curr to new state */                                 \
            t += h; /* advance time */                                          \
            for (int j = 0; j < VEC_SIZE; j++)                                  \
                y_curr[j] = sol_higher_ord[j];                                  \
                                                                                \
            /* Update solution structure */                                     \
            result->t[step] = t;                                                \
            for (int j = 0; j < VEC_SIZE; j++)                                  \
                /* write (separated for vectorization) */                       \
                result->y[step][j] = sol_higher_ord[j];                         \
        }                                                                       \
        /*  Otherwise, retry step */                                            \
                                                                                \
        /* Step size adjustment */                                              \
        /* Determine new step size (based on nth order local error) */          \
        h = 0.9 * h * pow(tol / err, 1.0 / ((double)SOLN_ORD));                 \
                                                                                \
        /* adjust timestep if needed to hit tf */                               \
        if (t + h > tf)                                                         \
            h = tf - t;                                                         \
    }                                                                           \
    /* record # of steps after finish */                                        \
    result->n = step;                                                           \
                                                                                \
    return result;

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "solver_weights.h"

typedef struct RKSolution
{
    double *t;             // timesteps
    double (*y)[VEC_SIZE]; // function value
    int n;                 // number of steps taken
} RKSolution;

RKSolution *rk89(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE]), double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK89_STAGES 16
#define RK89_ORD 9
    RK_METHOD(rk89_a, rk89_b, rk89_bh, rk89_c, RK89_ORD, RK89_STAGES)}

RKSolution *rk56(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE]), double t0, double tf, double y0[VEC_SIZE], double h0, double tol)
{
#define RK56_STAGES 9
#define RK56_ORD 6
    RK_METHOD(rk56_a, rk56_b, rk56_bh, rk56_c, 6.0, 9)
}

#endif
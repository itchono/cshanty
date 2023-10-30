// Numerical Runge Kutta Solvers

#ifndef SOLVERS_H
#define SOLVERS_H

#define VEC_SIZE 6 // dimension of vector (kept static at 6 for orbit prop)

// Big macro for defining general RK method
#define RK_METHOD_ADAPTIVE(weights_a, weights_b, weights_bh,                                  \
                           weights_c, SOLN_ORD, NUM_STAGES)                                   \
    if (tf <= t0)                                                                             \
        return NULL; /* error case */                                                         \
                                                                                              \
    /* Initial estimate for n = number of rows in solution */                                 \
    int n = (int)((tf - t0) / h0);                                                            \
    int step = 0;                                                                             \
                                                                                              \
    /* Allocate result struct */                                                              \
    RKSolution *result = (RKSolution *)malloc(sizeof(RKSolution));                            \
    result->y = (double(*)[VEC_SIZE])malloc(n * sizeof(double[VEC_SIZE]));                    \
    result->t = (double(*))malloc(n * sizeof(double));                                        \
    result->n_fev = 0;                                                                        \
    result->n_step_fail = 0;                                                                  \
                                                                                              \
    /* Prepare Step variables */                                                              \
    double h = h0;                  /* initial stepsize */                                    \
    double t = t0;                  /* intial time */                                         \
    double k[NUM_STAGES][VEC_SIZE]; /* k values; derivatives at each stage */                 \
    double y_stage[VEC_SIZE];       /* used to evaluate derivs at each stage */               \
    double y_curr[VEC_SIZE];        /* solution at beginning of current step */               \
                                                                                              \
    /* Init y0, t0 (step 0) */                                                                \
    for (int j = 0; j < VEC_SIZE; j++)                                                        \
    {                                                                                         \
        result->y[0][j] = y0[j];                                                              \
        y_curr[j] = y0[j];                                                                    \
    }                                                                                         \
    result->t[0] = t;                                                                         \
                                                                                              \
    /* Create Step Variables */                                                               \
                                                                                              \
    double err; /* magnitude of error, used for error estimate */                             \
    bool halt;                                                                                \
                                                                                              \
    /* zero-out k */                                                                          \
    for (int stage = 0; stage < NUM_STAGES; stage++)                                          \
        for (int j = 0; j < VEC_SIZE; j++)                                                    \
            k[stage][j] = 0;                                                                  \
                                                                                              \
    /* INTEGRATION STEPS */                                                                   \
    while ((t - h <= tf) && !halt)                                                            \
    {                                                                                         \
        /* Evaluate derivative at stage 0 */                                                  \
        f(t, y_curr, k[0], &halt);                                                            \
        result->n_fev++;                                                                      \
                                                                                              \
        /*Evaluate derivatives at stages 1-15 */                                              \
        for (int stage = 1; stage < NUM_STAGES; stage++)                                      \
        {                                                                                     \
            /* Prepare input vector */                                                        \
            for (int j = 0; j < VEC_SIZE; j++)                                                \
            {                                                                                 \
                y_stage[j] = y_curr[j]; /* take current sol */                                \
                for (int w = 0; w < stage; w++)                                               \
                    /* Add previous steps */                                                  \
                    y_stage[j] += h * k[w][j] * weights_a[stage][w];                          \
            }                                                                                 \
            /* evaluate next k */                                                             \
            f(t + h * weights_c[stage], y_stage, k[stage], &halt);                            \
            result->n_fev++;                                                                  \
        }                                                                                     \
                                                                                              \
        /* Calculate solution estimates */                                                    \
        err = 0;                                                                              \
        double sol_lower_ord[VEC_SIZE];                                                       \
        double sol_higher_ord[VEC_SIZE];                                                      \
                                                                                              \
        for (int j = 0; j < VEC_SIZE; j++)                                                    \
        {                                                                                     \
            sol_higher_ord[j] = y_curr[j];                                                    \
            sol_lower_ord[j] = y_curr[j];                                                     \
            for (int stage = 0; stage < NUM_STAGES; stage++)                                  \
            {                                                                                 \
                /* use b and b_h weights to calculate error */                                \
                sol_higher_ord[j] += h * weights_b[stage] * k[stage][j];                      \
                sol_lower_ord[j] += h * weights_bh[stage] * k[stage][j];                      \
            }                                                                                 \
            err += pow(sol_higher_ord[j] - sol_lower_ord[j], 2);                              \
        }                                                                                     \
        double old_err = err;                                                                 \
        err = sqrt(err);                                                                      \
        if (err != err)                                                                       \
        {                                                                                     \
            err = 0.9 * tol;                                                                  \
        }                                                                                     \
                                                                                              \
        if (err < tol)                                                                        \
        {                                                                                     \
            /* Accept step */                                                                 \
                                                                                              \
            /* Check array size and increase size if needed */                                \
            if ((step++) == n - 1)                                                            \
            {                                                                                 \
                n *= 2; /* double size of array and reallocate memory */                      \
                result->y = (double(*)[VEC_SIZE])realloc(                                     \
                    result->y, n * sizeof(double[VEC_SIZE]));                                 \
                result->t = (double(*))realloc(result->t, n * sizeof(double));                \
            }                                                                                 \
                                                                                              \
            /* Update t, y_curr to new state */                                               \
            t += h; /* advance time */                                                        \
            for (int j = 0; j < VEC_SIZE; j++)                                                \
                y_curr[j] = sol_higher_ord[j];                                                \
                                                                                              \
            /* Update solution structure */                                                   \
            result->t[step] = t;                                                              \
            for (int j = 0; j < VEC_SIZE; j++)                                                \
                /* write (separated for vectorization) */                                     \
                result->y[step][j] = sol_higher_ord[j];                                       \
        }                                                                                     \
        /*  Otherwise, retry step */                                                          \
        else                                                                                  \
        {                                                                                     \
            result->n_step_fail++;                                                            \
        }                                                                                     \
                                                                                              \
        /* Step size adjustment */                                                            \
        /* Determine new step size (based on nth order local error) */                        \
        h = 0.9 * h * pow(tol / err, 1.0 / ((double)SOLN_ORD));                               \
                                                                                              \
        /* adjust timestep if needed to hit tf */                                             \
        if (t + h > tf)                                                                       \
            h = tf - t;                                                                       \
    }                                                                                         \
    /* record # of steps after finish */                                                      \
    result->n = step;                                                                         \
    if (halt)                                                                                 \
    {                                                                                         \
        printf("CONVERGED!\n");                                                               \
    }                                                                                         \
    else                                                                                      \
    {                                                                                         \
        printf("NOT CONVERGED!, Time is %.4e, vs final which is %.4e, h = %.4e\n", t, tf, h); \
    }                                                                                         \
                                                                                              \
    return result;

#define RK_METHOD_FIXED_STEP(weights_a, weights_b, weights_bh,                                \
                             weights_c, SOLN_ORD, NUM_STAGES)                                 \
    if (tf <= t0)                                                                             \
        return NULL; /* error case */                                                         \
                                                                                              \
    /* Initial estimate for n = number of rows in solution */                                 \
    int n = (int)((tf - t0) / h);                                                             \
    int step = 0;                                                                             \
                                                                                              \
    /* Allocate result struct */                                                              \
    RKSolution *result = (RKSolution *)malloc(sizeof(RKSolution));                            \
    result->y = (double(*)[VEC_SIZE])malloc((2 * n) * sizeof(double[VEC_SIZE]));              \
    result->t = (double(*))malloc((2 * n) * sizeof(double));                                  \
    result->n_fev = 0;                                                                        \
    result->n_step_fail = 0;                                                                  \
                                                                                              \
    /* Prepare Step variables */                                                              \
    double t = t0;                  /* intial time */                                         \
    double k[NUM_STAGES][VEC_SIZE]; /* k values; derivatives at each stage */                 \
    double y_stage[VEC_SIZE];       /* used to evaluate derivs at each stage */               \
    double y_curr[VEC_SIZE];        /* solution at beginning of current step */               \
                                                                                              \
    /* Init y0, t0 (step 0) */                                                                \
    for (int j = 0; j < VEC_SIZE; j++)                                                        \
    {                                                                                         \
        result->y[0][j] = y0[j];                                                              \
        y_curr[j] = y0[j];                                                                    \
    }                                                                                         \
    result->t[0] = t;                                                                         \
                                                                                              \
    /* Create Step Variables */                                                               \
                                                                                              \
    bool halt;                                                                                \
                                                                                              \
    /* zero-out k */                                                                          \
    for (int stage = 0; stage < NUM_STAGES; stage++)                                          \
        for (int j = 0; j < VEC_SIZE; j++)                                                    \
            k[stage][j] = 0;                                                                  \
                                                                                              \
    /* INTEGRATION STEPS */                                                                   \
    while ((t - h <= tf) && !halt)                                                            \
    {                                                                                         \
        /* Evaluate derivative at stage 0 */                                                  \
        f(t, y_curr, k[0], &halt);                                                            \
        result->n_fev++;                                                                      \
                                                                                              \
        /*Evaluate derivatives at stages 1-15 */                                              \
        for (int stage = 1; stage < NUM_STAGES; stage++)                                      \
        {                                                                                     \
            /* Prepare input vector */                                                        \
            for (int j = 0; j < VEC_SIZE; j++)                                                \
            {                                                                                 \
                y_stage[j] = y_curr[j]; /* take current sol */                                \
                for (int w = 0; w < stage; w++)                                               \
                    /* Add previous steps */                                                  \
                    y_stage[j] += h * k[w][j] * weights_a[stage][w];                          \
            }                                                                                 \
            /* evaluate next k */                                                             \
            f(t + h * weights_c[stage], y_stage, k[stage], &halt);                            \
            result->n_fev++;                                                                  \
        }                                                                                     \
                                                                                              \
        /* Calculate solution estimates */                                                    \
        double sol_lower_ord[VEC_SIZE];                                                       \
        double sol_higher_ord[VEC_SIZE];                                                      \
                                                                                              \
        for (int j = 0; j < VEC_SIZE; j++)                                                    \
        {                                                                                     \
            sol_higher_ord[j] = y_curr[j];                                                    \
            sol_lower_ord[j] = y_curr[j];                                                     \
            for (int stage = 0; stage < NUM_STAGES; stage++)                                  \
            {                                                                                 \
                /* use b and b_h weights to calculate error */                                \
                sol_higher_ord[j] += h * weights_b[stage] * k[stage][j];                      \
                sol_lower_ord[j] += h * weights_bh[stage] * k[stage][j];                      \
            }                                                                                 \
        }                                                                                     \
                                                                                              \
        /* Update t, y_curr to new state */                                                   \
        step++;                                                                               \
        t += h; /* advance time */                                                            \
        for (int j = 0; j < VEC_SIZE; j++)                                                    \
            y_curr[j] = sol_higher_ord[j];                                                    \
                                                                                              \
        /* Update solution structure */                                                       \
        result->t[step] = t;                                                                  \
        for (int j = 0; j < VEC_SIZE; j++)                                                    \
            /* write (separated for vectorization) */                                         \
            result->y[step][j] = sol_higher_ord[j];                                           \
    }                                                                                         \
    /* record # of steps after finish */                                                      \
    result->n = step;                                                                         \
    if (halt)                                                                                 \
    {                                                                                         \
        printf("CONVERGED!\n");                                                               \
    }                                                                                         \
    else                                                                                      \
    {                                                                                         \
        printf("NOT CONVERGED!, Time is %.4e, vs final which is %.4e, h = %.4e\n", t, tf, h); \
    }                                                                                         \
                                                                                              \
    return result;

#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include "solver_weights.h"

typedef struct RKSolution
{
    double *t;             // timesteps
    double (*y)[VEC_SIZE]; // function value
    int n;                 // number of steps taken
    int n_fev;             // number of function evaluations
    int n_step_fail;       // number of failed steps
} RKSolution;

typedef RKSolution(*(*AdaptiveSolver)(void(double, double *, double *, bool *), double, double, double *, double, double));
typedef RKSolution(*(*FixedSolver)(void(double, double *, double *, bool *), double, double, double *, double));

// Professor Jim Verner's methods (rk89, rk78, rk67, rk56)

// robust 9th order RK method
RKSolution *rk89(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                 double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK89_STAGES 16
#define RK89_ORDER 9
    RK_METHOD_ADAPTIVE(rk89_a, rk89_b, rk89_bh, rk89_c, RK89_ORDER, RK89_STAGES)}

RKSolution *rk89_fixed(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                       double t0, double tf, double y0[VEC_SIZE], double h){
#define RK89_STAGES 16
#define RK89_ORDER 9
    RK_METHOD_FIXED_STEP(rk89_a, rk89_b, rk89_bh, rk89_c, RK89_ORDER, RK89_STAGES)}

// robust 8th order RK method
RKSolution *rk78(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                 double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK78_STAGES 13
#define RK78_ORDER 8
    RK_METHOD_ADAPTIVE(rk78_a, rk78_b, rk78_bh, rk78_c, RK78_ORDER, RK78_STAGES)}

// robust 7th order RK method
RKSolution *rk67(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                 double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK67_STAGES 10
#define RK67_ORDER 7
    RK_METHOD_ADAPTIVE(rk67_a, rk67_b, rk67_bh, rk67_c, RK67_ORDER, RK67_STAGES)}

// robust 6th order RK method
RKSolution *rk56(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                 double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK56_STAGES 9
#define RK56_ORDER 6
    RK_METHOD_ADAPTIVE(rk56_a, rk56_b, rk56_bh, rk56_c, RK56_ORDER, RK56_STAGES)}

// Professor Terry Feagin's methods (rk810, rk1012, rk1214)

// Less robust 10th order method
RKSolution *rk810(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                  double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK810_STAGES 17
#define RK810_ORDER 10
    RK_METHOD_ADAPTIVE(rk810_a, rk810_b, rk810_bh, rk810_c, RK810_ORDER, RK810_STAGES)}

// Less robust 12th order method
RKSolution *rk1012(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                   double t0, double tf, double y0[VEC_SIZE], double h0, double tol){
#define RK1012_STAGES 25
#define RK1012_ORDER 12
    RK_METHOD_ADAPTIVE(rk1012_a, rk1012_b, rk1012_bh, rk1012_c, RK1012_ORDER, RK1012_STAGES)}

// Less robust 14th order method
RKSolution *rk1214(void (*f)(double, double[VEC_SIZE], double[VEC_SIZE], bool *),
                   double t0, double tf, double y0[VEC_SIZE], double h0, double tol)
{
#define RK1214_STAGES 35
#define RK1214_ORDER 14
    RK_METHOD_ADAPTIVE(rk1214_a, rk1214_b, rk1214_bh, rk1214_c, RK1214_ORDER, RK1214_STAGES)
}

#endif
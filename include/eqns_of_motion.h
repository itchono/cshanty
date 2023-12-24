#ifndef EQNS_OF_MOTION_H
#define EQNS_OF_MOTION_H

#include <math.h>
#include "include/constants.h"

void gve_coeffs(double y[6], double A[6][3])
/**
 * @brief Calculate the coefficients of the Gauss variational equations
 *
 */
{
    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double q = 1 + f * cos(L) + g * sin(L);
    double leading_coeff = 1.0 / q * sqrt(p / mu);

    double pre_A[6][3] = {
        {0, 2 * p, 0},
        {q * sin(L),
         (q + 1) * cos(L) + f,
         -g * (h * sin(L) - k * cos(L))},
        {-q * cos(L),
         (q + 1) * sin(L) + g,
         f * (h * sin(L) - k * cos(L))},
        {0,
         0,
         cos(L) / 2 * (1 + h * h + k * k)},
        {0, 0, sin(L) / 2 * (1 + h * h + k * k)},
        {0,
         0,
         h * sin(L) - k * cos(L)}};

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A[i][j] = leading_coeff * pre_A[i][j];
        }
    }
}

void gauss_variational_eqns_mee(double t, double y[6], double dydt[6], double f_app[3])
{
    // Assume we get the unscaled value of p (i.e. in meters)

    double A[6][3];
    gve_coeffs(y, A);

    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double q = 1 + f * cos(L) + g * sin(L);

    double b[6] = {0,
                   0,
                   0,
                   0,
                   0,
                   q * q * sqrt(mu * p) / (p * p)};

    // y' = A*x + b
    for (int i = 0; i < 6; i++)
    {
        dydt[i] = b[i];

        for (int j = 0; j < 3; j++)
        {
            dydt[i] += A[i][j] * f_app[j];
        }
    }
}

#endif
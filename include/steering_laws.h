#ifndef STEERING_LAWS_H
#define STEERING_LAWS_H

#include <math.h>
#include "include/constants.h"
#include "include/rotations.h"
#include "include/math_utils.h"
#include "include/ephemeris.h"
#include "include/types.h"
#include "include/eqns_of_motion.h"

void approx_max_roc(double y[6], double maxroc[5])
/**
 * @brief Approximate maximum rates of change of each element across all steering angles and values of L.
 */
{
    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double q = 1 + f * cos(L) + g * sin(L);

    maxroc[0] = 2 * p / q * sqrt(p / mu);
    maxroc[1] = 2 * sqrt(p / mu); // approximations
    maxroc[2] = 2 * sqrt(p / mu); // approximations
    maxroc[3] = 1 / 2 * sqrt(p / mu) * (1 + h * h + k * k) / (sqrt(1 - g * g) + f);
    maxroc[4] = 1 / 2 * sqrt(p / mu) * (1 + h * h + k * k) / (sqrt(1 - f * f) + g);
}

void pe_penalty(double y[6], double pen_param, double rpmin, double *P, double dPdy[6])
{
    double p = y[0];
    double f = y[1];
    double g = y[2];

    double rp = p * (1 - sqrt(f * f + g * g)) / (1 - f * f - g * g);

    *P = exp(pen_param * (1 - rp / rpmin));

    dPdy[0] = (f * pen_param * p * exp(-pen_param * (p / (rpmin * (sqrt(f * f + g * g) + 1.0)) - 1.0)) * 1.0 / ((sqrt(f * f + g * g) + 1.0) * (sqrt(f * f + g * g) + 1.0)) * 1.0 / sqrt(f * f + g * g)) / rpmin;
    dPdy[1] = (g * pen_param * p * exp(-pen_param * (p / (rpmin * (sqrt(f * f + g * g) + 1.0)) - 1.0)) * 1.0 / ((sqrt(f * f + g * g) + 1.0) * (sqrt(f * f + g * g) + 1.0)) * 1.0 / sqrt(f * f + g * g)) / rpmin;
    dPdy[2] = -(pen_param * exp(-pen_param * (p / (rpmin * (sqrt(f * f + g * g) + 1.0)) - 1.0))) / (rpmin * (sqrt(f * f + g * g) + 1.0));
    dPdy[3] = 0.0;
    dPdy[4] = 0.0;
    dPdy[5] = 0.0;
}

void lyapunov_steering(double t, double y[6], ConfigStruct *cfg, double angles[2])
{
    double S[5] = {1. / 6378e3, 1, 1, 1, 1};
    double A[6][3];
    gve_coeffs(y, A);

    double d_oe_max[5];
    approx_max_roc(y, d_oe_max);

    double oe[5] = {y[0], y[1], y[2], y[3], y[4]};
    double oe_hat[5] = {cfg->y_target[0], cfg->y_target[1], cfg->y_target[2], cfg->y_target[3], cfg->y_target[4]};

    double P;
    double dPdoe[5];
    pe_penalty(y, cfg->penalty_param, cfg->min_pe, &P, dPdoe);

    double w_p = cfg->penalty_weight;

    double Xi[5]; // 5 x 1 vector
    for (int i = 0; i < 5; i++)
    {
        double Xi_penalty = dPdoe[i] * ((oe[i] - oe_hat[i]) / d_oe_max[i] * (oe[i] - oe_hat[i]) / d_oe_max[i]);
        double Xi_classic = 2 * (oe[i] - oe_hat[i]) / d_oe_max[i];

        Xi[i] = cfg->guidance_weights[i] * S[i] * (w_p * Xi_penalty + (1 + w_p * P) * Xi_classic);
    }

    double A_T[3][5]; // 3 x 5 matrix
    for (int i = 0; i < 5; i++)
    {
        // Take the first 5 rows of A and transpose them
        A_T[0][i] = A[i][0];
        A_T[1][i] = A[i][1];
        A_T[2][i] = A[i][2];
    }

    // final result is 3 x 1 vector achieved by A_T * Xi
    double d_Gamma_d_F[3];
    for (int i = 0; i < 3; i++)
    {
        d_Gamma_d_F[i] = 0;
        for (int j = 0; j < 5; j++)
        {
            d_Gamma_d_F[i] += A_T[i][j] * Xi[j];
        }
    }

    // GVE ordering is r, t, n; D1-3 is in order t, r, n
    double D1 = d_Gamma_d_F[1];
    double D2 = d_Gamma_d_F[0];
    double D3 = d_Gamma_d_F[2];

    angles[0] = atan2(-D2, -D1);
    angles[1] = atan2(-D3, sqrt(D1 * D1 + D2 * D2));
}

void ndf_heuristic(double t, double y[6], double ideal_angles[2], double adapted_angles[2])
{
    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double CIO[3][3];
    rot_inertial_LVLH(p, f, g, h, k, L, CIO);
    double COI[3][3];
    mat_transpose(CIO, COI);

    double lvlh[3];
    steering2lvlh(ideal_angles, lvlh);

    double n_star_i[3];
    mat_times_vec(CIO, lvlh, n_star_i);

    double sun_dir[3];
    sun_direction(t, sun_dir);

    double u_i[3] = {-sun_dir[0],
                     -sun_dir[1],
                     -sun_dir[2]};
    double c_cone_ang = vec_dot(n_star_i, u_i);

    // Perform adaptation

    // Sun-relative frame vectors
    double b_i[3];
    double n_i[3];
    vec_cross(n_star_i, u_i, n_i);
    vec_cross(u_i, n_i, b_i);

    // Target LVLH pointing vector
    double n_prime_o[3];

    if (c_cone_ang < cos(kappa_feathered))
    {
        // Feather the sail if we're below kappa_f, point in b_i
        mat_times_vec(COI, b_i, n_prime_o);
        lvlh2steering(n_prime_o, adapted_angles);
    }
    else if (c_cone_ang < cos(kappa_degraded))
    {
        // Degraded guidance, if we're between kappa_d and kappa_f vector
        // in plane of u_i and n_star, but orthogonal to u_i

        double n_prime_i[] = {cos(kappa_degraded) * u_i[0] + sin(kappa_degraded) * b_i[0],
                              cos(kappa_degraded) * u_i[1] + sin(kappa_degraded) * b_i[1],
                              cos(kappa_degraded) * u_i[2] + sin(kappa_degraded) * b_i[2]};
        mat_times_vec(COI, n_prime_i, n_prime_o);
        lvlh2steering(n_prime_o, adapted_angles);
    }
    else
    {
        adapted_angles[0] = ideal_angles[0];
        adapted_angles[1] = ideal_angles[1];
    }
}

#endif
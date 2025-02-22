#ifndef STEERING_LAWS_H
#define STEERING_LAWS_H

#include <math.h>
#include "constants.h"
#include "rotations.h"
#include "math_utils.h"
#include "ephemeris.h"
#include "types.h"
#include "eqns_of_motion.h"

void approx_max_roc(double y[6], double maxroc[5], double accel)
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

    maxroc[0] = 2.0 * p / q * sqrt(p / mu) * accel;
    maxroc[1] = 2.0 * sqrt(p / mu) * accel; // approximations
    maxroc[2] = 2.0 * sqrt(p / mu) * accel; // approximations
    maxroc[3] = 1.0 / 2.0 * sqrt(p / mu) * (1 + h * h + k * k) / (sqrt(1 - g * g) + f) * accel;
    maxroc[4] = 1.0 / 2.0 * sqrt(p / mu) * (1 + h * h + k * k) / (sqrt(1 - f * f) + g) * accel;
}

void max_roc_partials(double y[6], double partials[5], double accel)
{
    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double q = 1 + f * cos(L) + g * sin(L);

    partials[0] = 3.0 / q * sqrt(p / mu) * accel;
    partials[1] = 0;
    partials[2] = 0;
    partials[3] = sqrt(p / mu) * h / (sqrt(1 - g * g) + f) * accel;
    partials[4] = sqrt(p / mu) * k / (sqrt(1 - f * f) + g) * accel;
}

void pe_penalty(double y[6], double pen_param, double rpmin, double *P, double dPdy[5])
{
    double p = y[0];
    double f = y[1];
    double g = y[2];

    double rp = p * (1.0 - sqrt(f * f + g * g)) / (1 - f * f - g * g);

    *P = exp(pen_param * (1 - rp / rpmin));

    dPdy[0] = -(pen_param * exp(-pen_param * (p / (rpmin * (sqrt(f * f + g * g) + 1.0)) - 1.0))) / (rpmin * (sqrt(f * f + g * g) + 1.0));
    dPdy[1] = (f * pen_param * p * exp(-pen_param * (p / (rpmin * (sqrt(f * f + g * g) + 1.0)) - 1.0)) * 1.0 / ((sqrt(f * f + g * g) + 1.0) * (sqrt(f * f + g * g) + 1.0)) * 1.0 / sqrt(f * f + g * g)) / rpmin;
    dPdy[2] = (g * pen_param * p * exp(-pen_param * (p / (rpmin * (sqrt(f * f + g * g) + 1.0)) - 1.0)) * 1.0 / ((sqrt(f * f + g * g) + 1.0) * (sqrt(f * f + g * g) + 1.0)) * 1.0 / sqrt(f * f + g * g)) / rpmin;
    dPdy[3] = 0.0;
    dPdy[4] = 0.0;
}

void lyapunov_steering(double t, double y[6], ConfigStruct *cfg, double angles[2])
{
    double A[6][3];
    gve_coeffs(y, A);

    double accel = 2.0 * sail_p * sail_eta / cfg->sail_sigma;

    double d_oe_max[5];
    double partial_doe[5];
    approx_max_roc(y, d_oe_max, accel);
    max_roc_partials(y, partial_doe, accel);

    double oe[5] = {y[0], y[1], y[2], y[3], y[4]};
    double oe_hat[5] = {cfg->y_target[0], cfg->y_target[1], cfg->y_target[2], cfg->y_target[3], cfg->y_target[4]};

    double P;
    double dPdoe[5];
    pe_penalty(y, cfg->penalty_param, cfg->min_pe, &P, dPdoe);

    double W_p = cfg->penalty_weight;

    double Xi[5]; // 5 x 1 vector
    for (int i = 0; i < 5; i++)
    {
        double Xi_P = dPdoe[i] * ((oe[i] - oe_hat[i]) / d_oe_max[i] * (oe[i] - oe_hat[i]) / d_oe_max[i]);
        double Xi_Q = (oe[i] - oe_hat[i]) / (d_oe_max[i] * d_oe_max[i]);
        double Xi_R = -(oe[i] - oe_hat[i]) / (d_oe_max[i] * d_oe_max[i] * d_oe_max[i]) * partial_doe[i];

        Xi[i] = W_p * Xi_P + 2 * (1 + W_p * P) * (Xi_Q + Xi_R);
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
    double D[3] = {0, 0, 0};
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            D[i] += A_T[i][j] * cfg->guidance_weights[j] * Xi[j];
        }
    }

    // NAN check
    if (D[0] != D[0] || D[1] != D[1] || D[2] != D[2])
    {
        printf("NaN in steering law\n");
        angles[0] = 0;
        angles[1] = 0;
    }
    else
    {
        angles[0] = atan2(-D[0], -D[1]);
        angles[1] = atan2(-D[2], sqrt(D[0] * D[0] + D[1] * D[1]));
    }
}

void ndf_heuristic(double t, double y[6], double ideal_angles[2], ConfigStruct *cfg, double adapted_angles[2])
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

    if (c_cone_ang < cos(cfg->kappa_feathered))
    {
        // Feather the sail if we're below kappa_f, point in b_i
        mat_times_vec(COI, b_i, n_prime_o);
        lvlh2steering(n_prime_o, adapted_angles);
    }
    else if (c_cone_ang < cos(cfg->kappa_degraded))
    {
        // Degraded guidance, if we're between kappa_d and kappa_f vector
        // in plane of u_i and n_star, but orthogonal to u_i

        double n_prime_i[] = {cos(cfg->kappa_degraded) * u_i[0] + sin(cfg->kappa_degraded) * b_i[0],
                              cos(cfg->kappa_degraded) * u_i[1] + sin(cfg->kappa_degraded) * b_i[1],
                              cos(cfg->kappa_degraded) * u_i[2] + sin(cfg->kappa_degraded) * b_i[2]};
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
#ifndef STEERING_LAWS_H
#define STEERING_LAWS_H

#include <math.h>
#include "include/constants.h"
#include "rotations.h"
#include "math_utils.h"
#include "ephemeris.h"

void lyapunov_steering(double y[6], double y_tgt[5], double angles[2])
{
    // 1. Unpack vectors

    // Assume we get a properly scaled version of p
    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double p_hat = y_tgt[0];
    double f_hat = y_tgt[1];
    double g_hat = y_tgt[2];
    double h_hat = y_tgt[3];
    double k_hat = y_tgt[4];

    // 2. Compute Di factors

    double q = 1 + f * cos(L) + g * sin(L);

    // Scale down p by Earth radius
    double D1 = 2 * (p - p_hat) / r_earth +
                (f - f_hat) * ((q + 1) / q * cos(L) + f / q) +
                (g - g_hat) * ((q + 1) / q * sin(L) + g / q);
    double D2 = (f - f_hat) * sin(L) -
                (g - g_hat) * cos(L);
    double D3 = -g / q * (f - f_hat) * (h * sin(L) - k * cos(L)) +
                f / q * (g - g_hat) * (h * sin(L) - k * cos(L)) +
                2 * (sqrt(1 - g * g) + f) / q * cos(L) * (h - h_hat) +
                2 * (sqrt(1 - f * f) + g) / q * sin(L) * (k - k_hat);

    // 3. Compute optimal steering angles
    // alpha
    angles[0] = atan2(-D2, -D1);
    // beta
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
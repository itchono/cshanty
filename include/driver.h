#ifndef DRIVER_H
#define DRIVER_H

#include <stdbool.h>
#include "steering_laws.h"
#include "constants.h"
#include "eqns_of_motion.h"
#include "thrust.h"
#include "types.h"

#include <stdio.h>

void slyga_ode(double t, double y_sc[6], double dydt[6], bool *halt, bool *fault, ConfigStruct *cfg)
/**
 * @brief primary driving ODE for the solar sail sim
 *
 * @param t time
 * @param y_sc state vector
 * @param dydt derivative of state vector (output)
 * @param halt whether to halt the integration (output)
 * @param fault whether an error has occurred (output)
 * @param cfg configuration struct
 */
{
    // scaling
    double y[6];
    double S_ode[6] = {1. / cfg->y_target[0], 1, 1, 1, 1, 1};
    for (int i = 0; i < 6; i++)
    {
        y[i] = y_sc[i] / S_ode[i];
    }

    if ((fabs(y[1]) > 1) || (fabs(y[2]) > 1))
    {
        // eccentricity is greater than 1
        printf("ECCENTRICITY > 1: y[1] = %f, y[2] = %f\n", y[1], y[2]);
        *fault = true | *fault;
    }

    double ideal_angles[2];
    lyapunov_steering(t, y, cfg, ideal_angles);

    double actual_angles[2];
    ndf_heuristic(t, y, ideal_angles, cfg, actual_angles);

    double accel_o[3];
    sail_thrust(t, y, actual_angles, accel_o);

    double accel_norm = vec_norm(accel_o); // do something with this later for delta-v

    gauss_variational_eqns_mee(t, y, dydt, accel_o);

    // scale dydt
    for (int i = 0; i < 6; i++)
    {
        dydt[i] *= S_ode[i];
    }

    if (y[0] < r_earth && y[0] == y[0])
    {
        printf("plunged into the Earth (%f)\n", y[0]);
        *fault = true | *fault;
    }
    else
    {
        // calculate steering loss
        double err = 0;
        double S_guidance[5] = {1. / r_earth, 1, 1, 1, 1};
        for (int i = 0; i < 5; i++)
        {
            err += S_guidance[i] * (y[i] - cfg->y_target[i]) * (y[i] - cfg->y_target[i]);
        }
        err = sqrt(err);

        *halt = (bool)(err < cfg->guidance_tol) | *halt;
    }
    if (accel_norm != accel_norm)
    {
        printf("accel_norm is NaN\n");
        *fault = true | *fault;
    }
}

RKSolution *run_mission(ConfigStruct *cfg)
{
    ODESolver solver = *(cfg->solver);

    // pre-scale y0 by y_target
    cfg->y0[0] /= cfg->y_target[0];

    RKSolution *sol = solver(slyga_ode, cfg->t_span[0], cfg->t_span[1], cfg);

    // post-scale y by y_target
    for (int i = 0; i < sol->n; i++)
    {
        sol->y[i][0] *= cfg->y_target[0];
    }

    return sol;
}

#endif
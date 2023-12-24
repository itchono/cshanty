#ifndef DRIVER_H
#define DRIVER_H

#include <stdbool.h>
#include "include/steering_laws.h"
#include "include/constants.h"
#include "include/eqns_of_motion.h"
#include "include/thrust.h"
#include "include/types.h"

#include <stdio.h>

void slyga_ode(double t, double y[6], double dydt[6], bool *halt, ConfigStruct *cfg)
{

    double ideal_angles[2];
    lyapunov_steering(t, y, cfg, ideal_angles);

    double actual_angles[2];
    ndf_heuristic(t, y, ideal_angles, actual_angles);

    double accel_o[3];
    sail_thrust(t, y, actual_angles, accel_o);

    double accel_norm = vec_norm(accel_o); // do something with this later for delta-v

    gauss_variational_eqns_mee(t, y, dydt, accel_o);

    if (y[0] < 1 && y[0] == y[0])
    {
        printf("plunged into the Earth (%f)\n", y[0]);
        *halt = true;
    }
    else
    {
        // calculate steering loss
        double err = 0;
        double S[5] = {1. / 6378e3, 1, 1, 1, 1};
        for (int i = 0; i < 5; i++)
        {
            err += S[i] * (y[i] - cfg->y_target[i]) * (y[i] - cfg->y_target[i]);
        }
        err = sqrt(err);

        *halt = (bool)(err < cfg->guidance_tol);
    }
}

RKSolution *run_mission(ConfigStruct *cfg)
{
    // For now, run ODE unscaled
    ODESolver solver = *(cfg->solver);
    return solver(slyga_ode, cfg->t_span[0], cfg->t_span[1], cfg);
}

#endif
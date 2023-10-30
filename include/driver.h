#ifndef DRIVER_H
#define DRIVER_H

#include <stdbool.h>
#include "include/steering_laws.h"
#include "include/constants.h"
#include "include/eqns_of_motion.h"
#include "include/thrust.h"

#include <stdio.h>

extern double y_target[5];
extern double convergence_tol;

double convergence_tol = 0;
double y_target[5] = {0, 0, 0, 0, 0};

void slyga_ode(double t, double y[6], double dydt[6], bool *halt)
{

    // scaling
    double y_scaled[6] = {y[0] * r_earth,
                          y[1],
                          y[2],
                          y[3],
                          y[4],
                          y[5]};
    double y_tgt_scaled[6] = {y_target[0] * r_earth,
                              y_target[1],
                              y_target[2],
                              y_target[3],
                              y_target[4],
                              y_target[5]};

    double ideal_angles[2];
    lyapunov_steering(y_scaled, y_tgt_scaled, ideal_angles);

    double actual_angles[2];
    ndf_heuristic(t, y_scaled, ideal_angles, actual_angles);

    double accel_o[3];
    sail_thrust(t, y_scaled, actual_angles, accel_o);

    double y_p_unscaled[6];
    gauss_variational_eqns_mee(t, y_scaled, y_p_unscaled, accel_o);

    dydt[0] = y_p_unscaled[0] / r_earth;
    dydt[1] = y_p_unscaled[1];
    dydt[2] = y_p_unscaled[2];
    dydt[3] = y_p_unscaled[3];
    dydt[4] = y_p_unscaled[4];
    dydt[5] = y_p_unscaled[5];

    if (y[0] < 1 && y[0] == y[0])
    {
        printf("plunged into the Earth (%f)\n", y[0]);
        *halt = true;
    }
    else
    {
        // calculate steering loss
        double err = 0;
        for (int i = 0; i < 5; i++)
        {
            err += (y[i] - y_target[i]) * (y[i] - y_target[i]);
        }
        err = sqrt(err);

        *halt = (bool)(err < convergence_tol);
    }
}

#endif
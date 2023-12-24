#ifndef THRUST_H
#define THRUST_H

#include "include/converters.h"
#include "include/rotations.h"
#include "include/ephemeris.h"

void sail_thrust(double t, double y[6], double angles[2], double acceleration[3])
{
    double p = y[0];
    double f = y[1];
    double g = y[2];
    double h = y[3];
    double k = y[4];
    double L = y[5];

    double r_spacecraft_i[6];
    mee2cartesian(p, f, g, h, k, L, r_spacecraft_i);

    double CIO[3][3];
    rot_inertial_LVLH(p, f, g, h, k, L, CIO);
    double COI[3][3];
    mat_transpose(CIO, COI);

    double sc_dir_o[3];
    steering2lvlh(angles, sc_dir_o);

    double sc_dir_i[3];
    mat_times_vec(CIO, sc_dir_o, sc_dir_i);

    double sun_dir_i[3];
    sun_direction(t, sun_dir_i);

    double c_cone_ang = -vec_dot(sun_dir_i, sc_dir_i);

    // TODO add eclipse

    double efficiency = 2.0 * sail_p * sail_eta / sail_sigma;

    double thrust_sign;
    if (c_cone_ang > 0)
    {
        thrust_sign = 1;
    }
    else
    {
        thrust_sign = -1;
    }

    double prefactor = efficiency * thrust_sign * c_cone_ang * c_cone_ang;

    for (int i = 0; i < 3; i++)
        acceleration[i] = prefactor * sc_dir_o[i];
}

#endif
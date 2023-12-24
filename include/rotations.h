#ifndef ROTATIONS_H
#define ROTATIONS_H

#include "converters.h"
#include "math_utils.h"

void rot_inertial_LVLH(double p, double f, double g, double h, double k, double L, double CIO[3][3])
{
    double pos_vel[] = {0, 0, 0, 0, 0, 0};

    mee2cartesian(p, f, g, h, k, L, pos_vel);

    double mag_pos = vec_norm(pos_vel);
    double mag_vel = vec_norm(pos_vel + 3);

    double pos_u[3] = {pos_vel[0] / mag_pos,
                       pos_vel[1] / mag_pos,
                       pos_vel[2] / mag_pos};

    double vel_u[3] = {pos_vel[3] / mag_vel,
                       pos_vel[4] / mag_vel,
                       pos_vel[5] / mag_vel};

    double h_vec[] = {0, 0, 0};
    double b_vec[] = {0, 0, 0};
    vec_cross(pos_u, vel_u, h_vec);
    vec_cross(h_vec, pos_u, b_vec);

    for (int i = 0; i < 3; i++)
    {

        CIO[i][0] = pos_u[i];
        CIO[i][1] = b_vec[i];
        CIO[i][2] = h_vec[i];
    }
}

#endif
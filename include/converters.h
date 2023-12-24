#ifndef CONVERTERS_H
#define CONVERTERS_H

#include <math.h>
#include "include/constants.h"

void mee2cartesian(double p, double f, double g, double h, double k, double L, double pos_vel[6])
{
    // formulation from
    // https://spsweb.fltops.jpl.nasa.gov/portaldataops/mpg/MPG_Docs/Source%20Docs/EquinoctalElements-modified.pdf
    double alpha_sq = h * h - k * k;
    double s_sq = 1 + h * h + k * k;
    double q = 1 + f * cos(L) + g * sin(L);
    double r = p / q;

    double pos_fac = r / (s_sq);
    pos_vel[0] = pos_fac * (cos(L) + alpha_sq * cos(L) + 2 * h * k * sin(L));
    pos_vel[1] = pos_fac * (sin(L) - alpha_sq * sin(L) + 2 * h * k * cos(L));
    pos_vel[2] = pos_fac * (2 * (h * sin(L) - k * cos(L)));

    double vel_fac = 1.0 / (s_sq)*sqrt(mu / p);

    pos_vel[3] = vel_fac * (-(sin(L) + alpha_sq * sin(L) - 2 * h * k * cos(L) + g - 2 * f * h * k + alpha_sq * g));
    pos_vel[4] = vel_fac * (-(-cos(L) + alpha_sq * cos(L) + 2 * h * k * sin(L) - g + 2 * g * h * k + alpha_sq * f));
    pos_vel[5] = vel_fac * (2 * (h * cos(L) + k * sin(L) + f * h + g * k));
}

void steering2lvlh(double angles[2], double lvlh[3])
{
    double alpha = angles[0];
    double beta = angles[1];

    lvlh[0] = cos(beta) * sin(alpha);
    lvlh[1] = cos(beta) * cos(alpha);
    lvlh[2] = sin(beta);
}

void lvlh2steering(double lvlh[3], double angles[2])
{
    angles[0] = atan2(lvlh[0], lvlh[1]);
    angles[1] = atan2(lvlh[2], sqrt(lvlh[0] * lvlh[0] + lvlh[1] * lvlh[1]));
}
#endif
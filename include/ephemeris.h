#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include <math.h>
#include "constants.h"

double sun_angle(double t)
{
    const double t_sun = 31557600;
    return 2 * pi * t / t_sun;
}

void sun_direction(double t, double dir[3])
{
    double epsilon = 0.409087723; // obliquity of the ecliptic
    double lambda = sun_angle(t);
    dir[0] = cos(lambda);
    dir[1] = sin(lambda) * cos(epsilon);
    dir[2] = sin(lambda) * sin(epsilon);
}

#endif
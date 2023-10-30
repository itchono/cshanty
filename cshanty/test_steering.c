#include "include/steering_laws.h"
#include "include/driver.h"
#include <stdio.h>

int main()
{
    double y0[6] = {20000e3,
                    0.5,
                    -0.2,
                    0.5,
                    0,
                    0};

    double y_tgt[5] = {25000e3,
                       0.2,
                       0.5,
                       0,
                       0.3};

    double ideal_angles[2];

    lyapunov_steering(y0, y_tgt, ideal_angles);

    double actual_angles[2];
    ndf_heuristic(0, y0, ideal_angles, actual_angles);

    printf("alpha = %f, beta = %f\n", ideal_angles[0], ideal_angles[1]);
    printf("alpha = %f, beta = %f\n", actual_angles[0], actual_angles[1]);
}
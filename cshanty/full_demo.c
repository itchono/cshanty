#include <stdio.h>
#include "include/driver.h"
#include "include/constants.h"
#include "include/solvers.h"

int main()
{
    double y0[6] = {20000e3 / r_earth,
                    0.5,
                    -0.2,
                    0.5,
                    0,
                    0};

    y_target[0] = 25000e3 / r_earth;
    y_target[1] = 0.2;
    y_target[2] = 0.5;
    y_target[3] = 0;
    y_target[4] = 0.3;
    convergence_tol = 1e-3;

    AdaptiveSolver solver = rk67;

    RKSolution *sol = solver(slyga_ode, 0, 1e8, y0, 5, 1e-12);

    printf("n: %d\n", sol->n);
    printf("n_fev: %d\n", sol->n_fev);
    printf("n_step_fail: %d\n", sol->n_step_fail);
    printf("t_final: %f\n", sol->t[sol->n]);

    for (int i = 0; i < 6; i++)
    {
        if (i == 0)
        {
            sol->y[sol->n][i] *= r_earth;
        }

        printf("y[%d] = %f\n", i, sol->y[sol->n][i]);
    }

    free(sol->y);
    free(sol->t);
    free(sol);

    printf("Freed.\n");
}
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
    convergence_tol = 1e-2;

    FixedSolver solver = rk89_fixed;

    RKSolution *sol = solver(slyga_ode, 0, 1e8, y0, 1e2);

    printf("n: %d\n", sol->n);
    printf("n_fev: %d\n", sol->n_fev);
    printf("n_step_fail: %d\n", sol->n_step_fail);
    printf("t_final: %f\n", sol->t[sol->n]);
    printf("number of revolutions: %d\n", (int)(sol->y[sol->n][5] / (2 * pi)));

    FILE *fpt;

    fpt = fopen("MyFile.csv", "w+");

    for (int s = 0; s < sol->n; s++)
    {
        fprintf(fpt, "t = %.4e; {", sol->t[s]);
        for (int i = 0; i < 6; i++)
        {
            if (i == 0)
            {
                sol->y[s][i] *= r_earth;
            }

            fprintf(fpt, "%.4e, ", sol->y[s][i]);
        }
        fprintf(fpt, "}\n");
    }

    free(sol->y);
    free(sol->t);
    free(sol);

    printf("Freed.\n");
}
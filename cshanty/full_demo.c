#include <stdio.h>
#include "include/driver.h"
#include "include/constants.h"
#include "include/solvers.h"

int main()
{
    ConfigStruct cfg = {
        .y0 = {20000e3 / r_earth,
               0.5,
               -0.2,
               0.5,
               0,
               0},
        .y_target = {25000e3 / r_earth,
                     0.2,
                     0.5,
                     0,
                     0.3,
                     0},
        .propulsion_model = sail_thrust,
        .solver = rk89,
        .steering_law = lyapunov_steering,
        .t_span = {0, 1e8},
        .ode_rel_tol = 1e-4,
        .ode_h0 = 1e2,
        .guidance_tol = 1e-2,
        .guidance_weights = {1, 1, 1, 1, 1},
        .penalty_param = 1e-2,
        .min_pe = 0.1,
        .penalty_weight = 1e-2};

    RKSolution *sol = run_mission(&cfg);

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

    fclose(fpt);

    // free stuff
    free(sol->t);
    free(sol->y);
    free(sol);
    printf("done\n");
}

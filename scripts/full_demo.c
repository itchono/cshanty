#include <stdio.h>
#include "driver.h"
#include "constants.h"
#include "solvers.h"

int main()
{
    ConfigStruct cfg = {
        .y0 = {20000e3,
               0.5,
               -0.2,
               0.5,
               0,
               0},
        .y_target = {25000e3,
                     0.2,
                     0.5,
                     0,
                     0.3,
                     0},
        .propulsion_model = sail_thrust,
        .solver = rk89,
        .steering_law = lyapunov_steering,
        .t_span = {0, 1e8},
        .ode_rel_tol = 1e-6,
        .ode_h0 = 1e2,
        .guidance_tol = 5e-2,
        .guidance_weights = {1, 1, 1, 1, 1},
        .penalty_param = 1,
        .min_pe = 6878e3,
        .penalty_weight = 0,
        .kappa_degraded = pi / 180.0 * 70,
        .kappa_feathered = pi / 180.0 * 91};

    RKSolution *sol = run_mission(&cfg);

    printf("n: %d\n", sol->n);
    printf("n_fev: %d\n", sol->n_fev);
    printf("n_step_fail: %d\n", sol->n_step_fail);
    printf("t_final: %f\n", sol->t[sol->n - 1]);
    printf("number of revolutions: %d\n", (int)(sol->y[sol->n - 1][5] / (2 * pi)));

    FILE *fpt;

    fpt = fopen("MyFile.csv", "w+");

    for (int s = 0; s < sol->n; s++)
    {
        fprintf(fpt, "t = %.4e; {", sol->t[s]);
        for (int i = 0; i < 6; i++)
        {
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

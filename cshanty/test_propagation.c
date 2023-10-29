#include <stdio.h>

#include "include/eqns_of_motion.h"
#include "include/solvers.h"

void ode(double t, double *y, double *dydt)
{
    double f_app[3] = {0, 0.001, 0};
    gauss_variational_eqns_mee(t, y, dydt, f_app);
}

int main()
{

    double y0[6] = {20000e3,
                    0.5,
                    -0.2,
                    0.5,
                    0,
                    0};

    ODESolver solver = rk89;

    RKSolution *sol = solver(ode, 0, 1e8, y0, 5, 1e-5);

    printf("n: %d\n", sol->n);
    printf("n_fev: %d\n", sol->n_fev);
    printf("n_step_fail: %d\n", sol->n_step_fail);

    printf("t\t\ty\n");
    // for (int i = 0; i < sol->n; i++)
    // {
    //     printf("%.4e\t%.10e\n", sol->t[i], sol->y[i][0]);
    // }

    free(sol->y);
    free(sol->t);

    printf("Freed.\n");
}
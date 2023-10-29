#include "include/solvers.h"
#include <stdio.h>
#include <math.h>

void test_ode(double t, double *y, double *dydt)
{
    dydt[0] = -y[0];
}

double analytic_solution(double t)
{
    return exp(-t);
}

int main(int argc, char *argv[])
{
    double y0[] = {1, 0, 0, 0, 0, 0};

    if (argc == 2)
    {
        printf("The argument supplied is %s\n", argv[1]);
    }

    ODESolver solver = rk67;

    RKSolution *sol = solver(test_ode, 0, 10, y0, 5, 1e-5);
    printf("n: %d\n", sol->n);
    printf("n_fev: %d\n", sol->n_fev);
    printf("n_step_fail: %d\n", sol->n_step_fail);

    printf("t\t\ty\t\ty_analytic\terror\n");
    for (int i = 0; i < sol->n; i++)
    {
        printf("%.4e\t%.4e\t%.4e\t%.4e\n", sol->t[i], sol->y[i][0], analytic_solution(sol->t[i]), sol->y[i][0] - analytic_solution(sol->t[i]));
    }

    free(sol->y);
    free(sol->t);

    printf("Freed.\n");
}
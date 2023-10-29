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

int main()
{
    double y0[] = {1, 0, 0, 0, 0, 0};
    RKSolution *sol = rk56(test_ode, 0, 10, y0, 0.1, 1e-6);
    printf("n: %d\n", sol->n);

    printf("t\t\ty\t\ty_analytic\terror\n");
    for (int i = 0; i < sol->n; i++)
    {
        printf("%.4e\t%.4e\t%.4e\t%.4e\n", sol->t[i], sol->y[i][0], analytic_solution(sol->t[i]), sol->y[i][0] - analytic_solution(sol->t[i]));
    }
}
#include "include/solvers.h"
#include <stdio.h>
#include <math.h>

void test_ode(double t, double *y, double *dydt)
{
    dydt[0] = y[0];
}

double analytic_solution(double t)
{
    return exp(t);
}

int main()
{
    double y0[] = {1, 0, 0, 0, 0, 0};
    solution *sol = rk89(test_ode, 0, 1, y0, 0.1, 1);
    printf("n: %d\n", sol->n);

    printf("t\ty\ty_analytic\terror\n");
    for (int i = 0; i < sol->n; i++)
    {
        printf("%f\t%f\t%f\t%f\n", sol->t[i], sol->y[i][0], analytic_solution(sol->t[i]), sol->y[i][0] - analytic_solution(sol->t[i]));
    }
}
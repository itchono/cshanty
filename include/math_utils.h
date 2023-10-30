#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <math.h>

double vec_norm(double vec[3])
{
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

void vec_cross(double a[3], double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void mat_transpose(double A[3][3], double B[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            B[j][i] = A[i][j];
        }
    }
}

void mat_times_vec(double m[3][3], double a[3], double b[3])
{

    for (int i = 0; i < 3; i++)
    {
        b[i] = 0;
        for (int j = 0; j < 3; j++)
        {
            b[i] += m[i][j] * a[j];
        }
    }
}

double vec_dot(double a[3], double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

#endif
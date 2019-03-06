//
// Created by konstantin on 27.02.19.
//

#ifndef NEWFEM_LOCALMATRIX_H
#define NEWFEM_LOCALMATRIX_H

#include "math.h"

struct LocalMatrix
{
    void laplass(double ** K, double square, double * a, double * b)
    {
        K[0][0] = (b[0] * b[0] + a[0] * a[0]) / square / 4.0;
        K[0][1] = (b[1] * b[0] + a[1] * a[0]) / square / 4.0;
        K[0][2] = (b[2] * b[0] + a[2] * a[0]) / square / 4.0;

        K[1][0] = (b[0] * b[1] + a[0] * a[1]) / square / 4.0;
        K[1][1] = (b[1] * b[1] + a[1] * a[1]) / square / 4.0;
        K[1][2] = (b[2] * b[1] + a[2] * a[1]) / square / 4.0;

        K[2][0] = (b[0] * b[2] + a[0] * a[2]) / square / 4.0;
        K[2][1] = (b[1] * b[2] + a[1] * a[2]) / square / 4.0;
        K[2][2] = (b[2] * b[2] + a[2] * a[2]) / square / 4.0;
    }

    void mass(double ** K, double square)
    {
        K[0][0] = K[1][1] = K[2][2] = square / 6.;
        K[0][1] = K[0][2] = K[1][0] = K[1][2] = K[2][0] = K[2][1] = square / 12.;
    }

    void convectiveMembers(double ** K, double * f, double *g, double * a, double * b)
    {
        K[0][0] = b[0] * (f[0] / 12.) + b[0] * (f[1] / 24.) + b[0] * (f[2] / 24.)/**/ +(a[0] * g[0]) / 12. + (a[0] * g[1]) / 24. + (a[0] * g[2]) / 24.;
        K[0][1] = b[1] * (f[0] / 12.) + b[1] * (f[1] / 24.) + b[1] * (f[2] / 24.)/**/ +(a[1] * g[0]) / 12. + (a[1] * g[1]) / 24. + (a[1] * g[2]) / 24.;
        K[0][2] = b[2] * (f[0] / 12.) + b[2] * (f[1] / 24.) + b[2] * (f[2] / 24.)/**/ +(a[2] * g[0]) / 12. + (a[2] * g[1]) / 24. + (a[2] * g[2]) / 24.;

        K[1][0] = b[0] * (f[0] / 24.) + b[0] * (f[1] / 12.) + b[0] * (f[2] / 24.)/**/ +(a[0] * g[0]) / 24. + (a[0] * g[1]) / 12. + (a[0] * g[2]) / 24.;
        K[1][1] = b[1] * (f[0] / 24.) + b[1] * (f[1] / 12.) + b[1] * (f[2] / 24.)/**/ +(a[1] * g[0]) / 24. + (a[1] * g[1]) / 12. + (a[1] * g[2]) / 24.;
        K[1][2] = b[2] * (f[0] / 24.) + b[2] * (f[1] / 12.) + b[2] * (f[2] / 24.)/**/ +(a[2] * g[0]) / 24. + (a[2] * g[1]) / 12. + (a[2] * g[2]) / 24.;

        K[2][0] = b[0] * (f[0] / 24.) + b[0] * (f[1] / 24.) + b[0] * (f[2] / 12.)/**/ +(a[0] * g[0]) / 24. + (a[0] * g[1]) / 24. + (a[0] * g[2]) / 12.;
        K[2][1] = b[1] * (f[0] / 24.) + b[1] * (f[1] / 24.) + b[1] * (f[2] / 12.)/**/ +(a[1] * g[0]) / 24. + (a[1] * g[1]) / 24. + (a[1] * g[2]) / 12.;
        K[2][2] = b[2] * (f[0] / 24.) + b[2] * (f[1] / 24.) + b[2] * (f[2] / 12.)/**/ +(a[2] * g[0]) / 24. + (a[2] * g[1]) / 24. + (a[2] * g[2]) / 12.;
    }

    void derivative(double ** K, double *a_b)
    {
        K[0][0] = K[1][0] = K[2][0] = a_b[0] / 6.;
        K[0][1] = K[1][1] = K[2][1] = a_b[1] / 6.;
        K[0][2] = K[1][2] = K[2][2] = a_b[2] / 6.;
    }

    void convectiveMembersUpwind(double **K, double *a, double *b, double *u, double *v, double h, double k, double square)
    {
        if (h != 0)
        {
            for (int i = 0; i < 3; i++) // dc/dx
            {
                if (u[i] != 0)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        K[i][j] = (b[j] / 24.0) *
                                  ((2.0 * u[0] + u[1] + u[2]) +
                                   (u[i] / fabs(u[i])) *
                                   (((1.0 / tanh((fabs(u[i]) * h) / (2.0 * k)) - (2.0 * k) / (fabs(u[i]) * h)) * //alpha
                                     h * b[i]) / square) *
                                   (u[0] + u[1] + u[2]));
                    }
                } else
                {
                    for (int j = 0; j < 3; j++)
                    {
                        K[i][j] = (b[j] / 24.0) * (2.0 * u[0] + u[1] + u[2]);
                    }
                }
            }

            for (int i = 0; i < 3; i++)  // dc/dy
            {
                if (v[i] != 0)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        K[i][j] += (a[j] / 24.0) *
                                   ((2.0 * v[0] + v[1] + v[2]) +
                                    (v[i] / fabs(v[i])) *
                                    (((1.0 / tanh((fabs(v[i]) * h) / (2.0 * k)) - (2.0 * k) / (fabs(v[i]) * h)) * //alpha
                                      h * a[i]) / square) *
                                    (v[0] + v[1] + v[2]));
                    }
                } else
                {
                    for (int j = 0; j < 3; j++)
                    {
                        K[i][j] += (a[j] / 24.0) * (2.0 * v[0] + v[1] + v[2]);
                    }
                }
            }

        } else
        {
            LocalMatrix::convectiveMembers(K, u, v, a, b);
        }

    }
};

#endif //NEWFEM_LOCALMATRIX_H

//
// Created by konstantin on 27.02.19.
//

#ifndef NEWFEM_LOCALMATRIX_H
#define NEWFEM_LOCALMATRIX_H

#include "math.h"

struct Supg
{
    void laplass(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix)
    {
        double uAvg = (u[0] + u[1] + u[2]) / 3.0;
        double vAvg = (v[0] + v[1] + v[2]) / 3.0;
        double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

        if (fabs(uMod) > 0.0001)
        {
            double pe = (uMod * h) / (2.0 * k);
            double tau = (h / (2.0 * uMod));

            if (k != 0)
                tau *= (1.0 / tanh(pe) - 1.0 / pe);

            matrix[0][0] = ((b[0]*b[0])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[0] + vAvg * a[0]));
            matrix[0][1] = ((b[0]*b[1])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[0] + vAvg * a[0]));
            matrix[0][2] = ((b[0]*b[2])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[0] + vAvg * a[0]));

            matrix[1][0] = ((b[1]*b[0])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[1] + vAvg * a[1]));
            matrix[1][1] = ((b[1]*b[1])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[1] + vAvg * a[1]));
            matrix[1][2] = ((b[1]*b[2])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[1] + vAvg * a[1]));

            matrix[2][0] = ((b[2]*b[0])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[2] + vAvg * a[2]));
            matrix[2][1] = ((b[2]*b[1])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[2] + vAvg * a[2]));
            matrix[2][2] = ((b[2]*b[2])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[2] + vAvg * a[2]));
            //////////////////////////////////////////////////////////////////////////////////////////////////////////
            matrix[0][0] += ((a[0]*a[0])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[0] + vAvg * a[0]));
            matrix[0][1] += ((a[0]*a[1])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[0] + vAvg * a[0]));
            matrix[0][2] += ((a[0]*a[2])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[0] + vAvg * a[0]));

            matrix[1][0] += ((a[1]*a[0])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[1] + vAvg * a[1]));
            matrix[1][1] += ((a[1]*a[1])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[1] + vAvg * a[1]));
            matrix[1][2] += ((a[1]*a[2])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[1] + vAvg * a[1]));

            matrix[2][0] += ((a[2]*a[0])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[2] + vAvg * a[2]));
            matrix[2][1] += ((a[2]*a[1])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[2] + vAvg * a[2]));
            matrix[2][2] += ((a[2]*a[2])/(4.0 * square)) * (1.0 / 3.0 + (tau / (2.0 * square)) * (uAvg * b[2] + vAvg * a[2]));
        } else
        {
            matrix[0][0] = (b[0] * b[0] + a[0] * a[0]) / square / 4.0;
            matrix[0][1] = (b[1] * b[0] + a[1] * a[0]) / square / 4.0;
            matrix[0][2] = (b[2] * b[0] + a[2] * a[0]) / square / 4.0;

            matrix[1][0] = (b[0] * b[1] + a[0] * a[1]) / square / 4.0;
            matrix[1][1] = (b[1] * b[1] + a[1] * a[1]) / square / 4.0;
            matrix[1][2] = (b[2] * b[1] + a[2] * a[1]) / square / 4.0;

            matrix[2][0] = (b[0] * b[2] + a[0] * a[2]) / square / 4.0;
            matrix[2][1] = (b[1] * b[2] + a[1] * a[2]) / square / 4.0;
            matrix[2][2] = (b[2] * b[2] + a[2] * a[2]) / square / 4.0;
        }
    }

    void convectiveMembers(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix)
    {
        double uAvg = (u[0] + u[1] + u[2]) / 3.0;
        double vAvg = (v[0] + v[1] + v[2]) / 3.0;
        double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

        if (fabs(uMod) > 0.0001) //1!
        {
            double alpha = 1;
            if (k != 0)
                alpha = (1.0 / tanh((uMod * h) / (2 * k))) - (2.0 * k) / (uMod * h);

            matrix[0][0] = b[0] * (2.0*u[0] + u [1] + u[2]) / 24.0 + uAvg * alpha * h * b[0] * (uAvg * b[0] + vAvg * a[0]) / (8.0 * square * uMod);
            matrix[0][1] = b[1] * (2.0*u[0] + u [1] + u[2]) / 24.0 + uAvg * alpha * h * b[1] * (uAvg * b[0] + vAvg * a[0]) / (8.0 * square * uMod);
            matrix[0][2] = b[2] * (2.0*u[0] + u [1] + u[2]) / 24.0 + uAvg * alpha * h * b[2] * (uAvg * b[0] + vAvg * a[0]) / (8.0 * square * uMod);

            matrix[1][0] = b[0] * (u[0] + 2.0*u [1] + u[2]) / 24.0 + uAvg * alpha * h * b[0] * (uAvg * b[1] + vAvg * a[1]) / (8.0 * square * uMod);
            matrix[1][1] = b[1] * (u[0] + 2.0*u [1] + u[2]) / 24.0 + uAvg * alpha * h * b[1] * (uAvg * b[1] + vAvg * a[1]) / (8.0 * square * uMod);
            matrix[1][2] = b[2] * (u[0] + 2.0*u [1] + u[2]) / 24.0 + uAvg * alpha * h * b[2] * (uAvg * b[1] + vAvg * a[1]) / (8.0 * square * uMod);

            matrix[2][0] = b[0] * (u[0] + u [1] + 2.0*u[2]) / 24.0 + uAvg * alpha * h * b[0] * (uAvg * b[2] + vAvg * a[2]) / (8.0 * square * uMod);
            matrix[2][1] = b[1] * (u[0] + u [1] + 2.0*u[2]) / 24.0 + uAvg * alpha * h * b[1] * (uAvg * b[2] + vAvg * a[2]) / (8.0 * square * uMod);
            matrix[2][2] = b[2] * (u[0] + u [1] + 2.0*u[2]) / 24.0 + uAvg * alpha * h * b[2] * (uAvg * b[2] + vAvg * a[2]) / (8.0 * square * uMod);



            matrix[0][0] += a[0] * (2.0*v[0] + v[1] + v[2]) / 24.0 + vAvg * alpha * h * a[0] * (uAvg * b[0] + vAvg * a[0]) / (8.0 * square * uMod);
            matrix[0][1] += a[1] * (2.0*v[0] + v[1] + v[2]) / 24.0 + vAvg * alpha * h * a[1] * (uAvg * b[0] + vAvg * a[0]) / (8.0 * square * uMod);
            matrix[0][2] += a[2] * (2.0*v[0] + v[1] + v[2]) / 24.0 + vAvg * alpha * h * a[2] * (uAvg * b[0] + vAvg * a[0]) / (8.0 * square * uMod);

            matrix[1][0] += a[0] * (v[0] + 2.0*v[1] + v[2]) / 24.0 + vAvg * alpha * h * a[0] * (uAvg * b[1] + vAvg * a[1]) / (8.0 * square * uMod);
            matrix[1][1] += a[1] * (v[0] + 2.0*v[1] + v[2]) / 24.0 + vAvg * alpha * h * a[1] * (uAvg * b[1] + vAvg * a[1]) / (8.0 * square * uMod);
            matrix[1][2] += a[2] * (v[0] + 2.0*v[1] + v[2]) / 24.0 + vAvg * alpha * h * a[2] * (uAvg * b[1] + vAvg * a[1]) / (8.0 * square * uMod);

            matrix[2][0] += a[0] * (v[0] + v[1] + 2.0*v[2]) / 24.0 + vAvg * alpha * h * a[0] * (uAvg * b[2] + vAvg * a[2]) / (8.0 * square * uMod);
            matrix[2][1] += a[1] * (v[0] + v[1] + 2.0*v[2]) / 24.0 + vAvg * alpha * h * a[1] * (uAvg * b[2] + vAvg * a[2]) / (8.0 * square * uMod);
            matrix[2][2] += a[2] * (v[0] + v[1] + 2.0*v[2]) / 24.0 + vAvg * alpha * h * a[2] * (uAvg * b[2] + vAvg * a[2]) / (8.0 * square * uMod);
        } else
        {
            matrix[0][0] = 0;
            matrix[0][1] = 0;
            matrix[0][2] = 0;

            matrix[1][0] = 0;
            matrix[1][1] = 0;
            matrix[1][2] = 0;

            matrix[2][0] = 0;
            matrix[2][1] = 0;
            matrix[2][2] = 0;
        }
    }

    void mass(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix)
    {
        double uAvg = (u[0] + u[1] + u[2]) / 3.0;
        double vAvg = (v[0] + v[1] + v[2]) / 3.0;
        double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

        if (fabs(uMod) > 0.0001)
        {
            double alpha = 1;
            if (k != 0)
                alpha = (1.0 / tanh((uMod * h) / (2 * k))) - (2.0 * k) / (uMod * h);

            matrix[0][0] = square / 6.0  + (alpha * h * (uAvg * b[0] + vAvg * a[0])) / (12.0 * uMod);
            matrix[0][1] = square / 12.0 + (alpha * h * (uAvg * b[0] + vAvg * a[0])) / (12.0 * uMod);
            matrix[0][2] = square / 12.0 + (alpha * h * (uAvg * b[0] + vAvg * a[0])) / (12.0 * uMod);

            matrix[1][0] = square / 12.0 + (alpha * h * (uAvg * b[1] + vAvg * a[1])) / (12.0 * uMod);
            matrix[1][1] = square / 6.0  + (alpha * h * (uAvg * b[1] + vAvg * a[1])) / (12.0 * uMod);
            matrix[1][2] = square / 12.0 + (alpha * h * (uAvg * b[1] + vAvg * a[1])) / (12.0 * uMod);

            matrix[2][0] = square / 12.0 + (alpha * h * (uAvg * b[2] + vAvg * a[2])) / (12.0 * uMod);
            matrix[2][1] = square / 12.0 + (alpha * h * (uAvg * b[2] + vAvg * a[2])) / (12.0 * uMod);
            matrix[2][2] = square / 6.0  + (alpha * h * (uAvg * b[2] + vAvg * a[2])) / (12.0 * uMod);
        } else
        {
            matrix[0][0] = square / 6.0;//  + (alpha * h * (uAvg * b[0] + vAvg * a[0])) / (12.0 * uMod);
            matrix[0][1] = square / 12.0;// + (alpha * h * (uAvg * b[0] + vAvg * a[0])) / (12.0 * uMod);
            matrix[0][2] = square / 12.0;// + (alpha * h * (uAvg * b[0] + vAvg * a[0])) / (12.0 * uMod);

            matrix[1][0] = square / 12.0;// + (alpha * h * (uAvg * b[1] + vAvg * a[1])) / (12.0 * uMod);
            matrix[1][1] = square / 6.0;//  + (alpha * h * (uAvg * b[1] + vAvg * a[1])) / (12.0 * uMod);
            matrix[1][2] = square / 12.0;// + (alpha * h * (uAvg * b[1] + vAvg * a[1])) / (12.0 * uMod);

            matrix[2][0] = square / 12.0;// + (alpha * h * (uAvg * b[2] + vAvg * a[2])) / (12.0 * uMod);
            matrix[2][1] = square / 12.0;// + (alpha * h * (uAvg * b[2] + vAvg * a[2])) / (12.0 * uMod);
            matrix[2][2] = square / 6.0;//  + (alpha * h * (uAvg * b[2] + vAvg * a[2])) / (12.0 * uMod);
        }
    }
};

struct LocalMatrix
{

    Supg supg;

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

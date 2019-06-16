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

struct TwoPhase
{

    void derPPP(double *a, double *rho, double *p, double *vector)
    {
        double temp = (p[0] * a[0] + p[1] * a[1] + p[2] * a[2]) / 24.0;

        vector[0] = temp * (2.0 / rho[0] + 1.0 / rho[1] + 1.0 / rho[2]);
        vector[1] = temp * (1.0 / rho[0] + 2.0 / rho[1] + 1.0 / rho[2]);
        vector[2] = temp * (1.0 / rho[0] + 1.0 / rho[1] + 2.0 / rho[2]);
    }

    void massRho(double square, double *rho, double **matrix)
    {
        matrix[0][0] = square * (rho[0] / 10.0 + rho[1] / 30.0 + rho[2] / 30.0);
        matrix[0][1] = square * (rho[0] / 30.0 + rho[1] / 30.0 + rho[2] / 60.0);
        matrix[0][2] = square * (rho[0] / 30.0 + rho[1] / 60.0 + rho[2] / 30.0);

        matrix[1][0] = square * (rho[0] / 30.0 + rho[1] / 30.0 + rho[2] / 60.0);
        matrix[1][1] = square * (rho[0] / 30.0 + rho[1] / 10.0 + rho[2] / 30.0);
        matrix[1][2] = square * (rho[0] / 60.0 + rho[1] / 30.0 + rho[2] / 30.0);

        matrix[2][0] = square * (rho[0] / 30.0 + rho[1] / 60.0 + rho[2] / 30.0);
        matrix[2][1] = square * (rho[0] / 60.0 + rho[1] / 30.0 + rho[2] / 30.0);
        matrix[2][2] = square * (rho[0] / 30.0 + rho[1] / 30.0 + rho[2] / 10.0);
    }

    void convectiveMembersRho(double *a, double *b, double *u, double *v, double *rho, double square, double **matrix)
    {
        double uTmp0 = (rho[0]*(3.0 * u[0] + u[1] + u[2])+rho[1]*(u[0] + u[1] + u[2] / 2.0) + rho[2] * (u[0] + u[1] / 2.0 + u[2])) / 60.0;
        double uTmp1 = (rho[0]*(u[0] + u[1] + u[2] / 2.0)+rho[1]*(u[0] + 3.0 * u[1] + u[2]) + rho[2] * (u[0] / 2.0 + u[1] + u[2])) / 60.0;
        double uTmp2 = (rho[0]*(u[0] + u[1] / 2.0 + u[2])+rho[1]*(u[0] / 2.0 + u[1] + u[2]) + rho[2] * (u[0] + u[1] + 3.0 * u[2])) / 60.0;

        double vTmp0 = (rho[0]*(3.0 * v[0] + v[1] + v[2])+rho[1]*(v[0] + v[1] + v[2] / 2.0) + rho[2] * (v[0] + v[1] / 2.0 + v[2])) / 60.0;
        double vTmp1 = (rho[0]*(v[0] + v[1] + v[2] / 2.0)+rho[1]*(v[0] + 3.0 * v[1] + v[2]) + rho[2] * (v[0] / 2.0 + v[1] + v[2])) / 60.0;
        double vTmp2 = (rho[0]*(v[0] + v[1] / 2.0 + v[2])+rho[1]*(v[0] / 2.0 + v[1] + v[2]) + rho[2] * (v[0] + v[1] + 3.0 * v[2])) / 60.0;

        matrix[0][0] = b[0] * uTmp0 + a[0] * vTmp0;
        matrix[0][1] = b[1] * uTmp0 + a[1] * vTmp0;
        matrix[0][2] = b[2] * uTmp0 + a[2] * vTmp0;

        matrix[1][0] = b[0] * uTmp2 + a[0] * vTmp1;
        matrix[1][1] = b[1] * uTmp2 + a[1] * vTmp1;
        matrix[1][2] = b[2] * uTmp2 + a[2] * vTmp1;

        matrix[2][0] = b[0] * uTmp2 + a[0] * vTmp2;
        matrix[2][1] = b[1] * uTmp2 + a[1] * vTmp2;
        matrix[2][2] = b[2] * uTmp2 + a[2] * vTmp2;
    }

    void laplassP(double *a, double *b, double *rho, double square, double **matrix)
    {
        double temp = (1.0 / rho[0] + 1.0 / rho[1] + 1.0 / rho[2]) / (12.0 * square);

        matrix[0][0] = (b[0] * b[0] + a[0]* a[0]) * temp;
        matrix[0][1] = (b[0] * b[1] + a[0]* a[1]) * temp;
        matrix[0][2] = (b[0] * b[2] + a[0]* a[2]) * temp;

        matrix[1][0] = (b[1] * b[0] + a[1]* a[0]) * temp;
        matrix[1][1] = (b[1] * b[1] + a[1]* a[1]) * temp;
        matrix[1][2] = (b[1] * b[2] + a[1]* a[2]) * temp;

        matrix[2][0] = (b[2] * b[0] + a[2]* a[0]) * temp;
        matrix[2][1] = (b[2] * b[1] + a[2]* a[1]) * temp;
        matrix[2][2] = (b[2] * b[2] + a[2]* a[2]) * temp;

//        double temp = 1.0 / (4.0 * square);
//
//        matrix[0][0] = temp * (b[0] * b[0] + a[0]* a[0]) / rho[0];
//        matrix[0][1] = temp * (b[0] * b[1] + a[0]* a[1]) / rho[1];
//        matrix[0][2] = temp * (b[0] * b[2] + a[0]* a[2]) / rho[2];
//
//        matrix[1][0] = temp * (b[1] * b[0] + a[1]* a[0]) / rho[0];
//        matrix[1][1] = temp * (b[1] * b[1] + a[1]* a[1]) / rho[1];
//        matrix[1][2] = temp * (b[1] * b[2] + a[1]* a[2]) / rho[2];
//
//        matrix[2][0] = temp * (b[2] * b[0] + a[2]* a[0]) / rho[0];
//        matrix[2][1] = temp * (b[2] * b[1] + a[2]* a[1]) / rho[1];
//        matrix[2][2] = temp * (b[2] * b[2] + a[2]* a[2]) / rho[2];

//        double temp = (b[0] / rho[0] + b[1] / rho[1] + b[2] / rho[2]) / (12.0 * square);
//
//        matrix[0][0] = b[0] * temp;
//        matrix[0][1] = b[1] * temp;
//        matrix[0][2] = b[2] * temp;
//
//        matrix[1][0] = b[0] * temp;
//        matrix[1][1] = b[1] * temp;
//        matrix[1][2] = b[2] * temp;
//
//        matrix[2][0] = b[0] * temp;
//        matrix[2][1] = b[1] * temp;
//        matrix[2][2] = b[2] * temp;
//
//        temp = (a[0] / rho[0] + a[1] / rho[1] + a[2] / rho[2]) / (12.0 * square);
//
//        matrix[0][0] += a[0] * temp;
//        matrix[0][1] += a[1] * temp;
//        matrix[0][2] += a[2] * temp;
//
//        matrix[1][0] += a[0] * temp;
//        matrix[1][1] += a[1] * temp;
//        matrix[1][2] += a[2] * temp;
//
//        matrix[2][0] += a[0] * temp;
//        matrix[2][1] += a[1] * temp;
//        matrix[2][2] += a[2] * temp;
//
//        temp = - 1.0 / (4.0 * square);
//        double Atemp = (b[0] / rho[0] + b[1] / rho[1] + b[2] / rho[2]) / 3.0;
//
//        matrix[0][0] += temp * (b[0] * (Atemp + b[0] / rho[0]));
//        matrix[0][1] += temp * (b[1] * (Atemp + b[0] / rho[0]));
//        matrix[0][2] += temp * (b[2] * (Atemp + b[0] / rho[0]));
//
//        matrix[1][0] += temp * (b[0] * (Atemp + b[1] / rho[1]));
//        matrix[1][1] += temp * (b[1] * (Atemp + b[1] / rho[1]));
//        matrix[1][2] += temp * (b[2] * (Atemp + b[1] / rho[1]));
//
//        matrix[2][0] += temp * (b[0] * (Atemp + b[2] / rho[2]));
//        matrix[2][1] += temp * (b[1] * (Atemp + b[2] / rho[2]));
//        matrix[2][2] += temp * (b[2] * (Atemp + b[2] / rho[2]));
//
//        Atemp = (a[0] / rho[0] + a[1] / rho[1] + a[2] / rho[2]) / 3.0;
//
//        matrix[0][0] += temp * (a[0] * (Atemp + a[0] / rho[0]));
//        matrix[0][1] += temp * (a[1] * (Atemp + a[0] / rho[0]));
//        matrix[0][2] += temp * (a[2] * (Atemp + a[0] / rho[0]));
//
//        matrix[1][0] += temp * (a[0] * (Atemp + a[1] / rho[1]));
//        matrix[1][1] += temp * (a[1] * (Atemp + a[1] / rho[1]));
//        matrix[1][2] += temp * (a[2] * (Atemp + a[1] / rho[1]));
//
//        matrix[2][0] += temp * (a[0] * (Atemp + a[2] / rho[2]));
//        matrix[2][1] += temp * (a[1] * (Atemp + a[2] / rho[2]));
//        matrix[2][2] += temp * (a[2] * (Atemp + a[2] / rho[2]));
//
//        matrix[0][0] *=-1;
//        matrix[0][1] *=-1;
//        matrix[0][2] *=-1;
//
//        matrix[1][0] *=-1;
//        matrix[1][1] *=-1;
//        matrix[1][2] *=-1;
//
//        matrix[2][0] *=-1;
//        matrix[2][1] *=-1;
//        matrix[2][2] *=-1;
    }

    void laplassMuX(double *a, double *b, double *mu, double *v, double square, double *vector, double **matrix)
    {
        double temp = (mu[0] + mu[1] + mu[2]) / (6.0 * square);

        matrix[0][0] = temp * b[0] * b[0];
        matrix[0][1] = temp * b[0] * b[1];
        matrix[0][2] = temp * b[0] * b[2];

        matrix[1][0] = temp * b[1] * b[0];
        matrix[1][1] = temp * b[1] * b[1];
        matrix[1][2] = temp * b[1] * b[2];

        matrix[2][0] = temp * b[2] * b[0];
        matrix[2][1] = temp * b[2] * b[1];
        matrix[2][2] = temp * b[2] * b[2];

        temp /= 2.0;

        matrix[0][0] += temp * a[0] * a[0];
        matrix[0][1] += temp * a[0] * a[1];
        matrix[0][2] += temp * a[0] * a[2];

        matrix[1][0] += temp * a[1] * a[0];
        matrix[1][1] += temp * a[1] * a[1];
        matrix[1][2] += temp * a[1] * a[2];

        matrix[2][0] += temp * a[2] * a[0];
        matrix[2][1] += temp * a[2] * a[1];
        matrix[2][2] += temp * a[2] * a[2];

        vector[0] = temp * a[0] * (b[0] * v[0] + b[1] * v[1] + b[2] * v[2]);
        vector[1] = temp * a[1] * (b[0] * v[0] + b[1] * v[1] + b[2] * v[2]);
        vector[2] = temp * a[2] * (b[0] * v[0] + b[1] * v[1] + b[2] * v[2]);
    }

    void laplassMuY(double *a, double *b, double *mu, double *u, double square, double *vector, double **matrix)
    {
        double temp = (mu[0] + mu[1] + mu[2]) / (6.0 * square);

        matrix[0][0] = temp * a[0] * a[0];
        matrix[0][1] = temp * a[0] * a[1];
        matrix[0][2] = temp * a[0] * a[2];

        matrix[1][0] = temp * a[1] * a[0];
        matrix[1][1] = temp * a[1] * a[1];
        matrix[1][2] = temp * a[1] * a[2];

        matrix[2][0] = temp * a[2] * a[0];
        matrix[2][1] = temp * a[2] * a[1];
        matrix[2][2] = temp * a[2] * a[2];

        temp /= 2.0;

        matrix[0][0] += temp * b[0] * b[0];
        matrix[0][1] += temp * b[0] * b[1];
        matrix[0][2] += temp * b[0] * b[2];

        matrix[1][0] += temp * b[1] * b[0];
        matrix[1][1] += temp * b[1] * b[1];
        matrix[1][2] += temp * b[1] * b[2];

        matrix[2][0] += temp * b[2] * b[0];
        matrix[2][1] += temp * b[2] * b[1];
        matrix[2][2] += temp * b[2] * b[2];

        vector[0] = temp * b[0] * (a[0] * u[0] + a[1] * u[1] + a[2] * u[2]);
        vector[1] = temp * b[1] * (a[0] * u[0] + a[1] * u[1] + a[2] * u[2]);
        vector[2] = temp * b[2] * (a[0] * u[0] + a[1] * u[1] + a[2] * u[2]);
    }
};

struct LocalMatrix
{

    Supg supg;
    TwoPhase twoPhase;

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

    void laplassX(double ** K, double square, double * a, double * b)
    {
        K[0][0] = (b[0] * b[0]) / square / 4.0;
        K[0][1] = (b[1] * b[0]) / square / 4.0;
        K[0][2] = (b[2] * b[0]) / square / 4.0;

        K[1][0] = (b[0] * b[1]) / square / 4.0;
        K[1][1] = (b[1] * b[1]) / square / 4.0;
        K[1][2] = (b[2] * b[1]) / square / 4.0;

        K[2][0] = (b[0] * b[2]) / square / 4.0;
        K[2][1] = (b[1] * b[2]) / square / 4.0;
        K[2][2] = (b[2] * b[2]) / square / 4.0;
    }

    void laplassY(double ** K, double square, double * a, double * b)
    {
        K[0][0] = (a[0] * a[0]) / square / 4.0;
        K[0][1] = (a[1] * a[0]) / square / 4.0;
        K[0][2] = (a[2] * a[0]) / square / 4.0;

        K[1][0] = (a[0] * a[1]) / square / 4.0;
        K[1][1] = (a[1] * a[1]) / square / 4.0;
        K[1][2] = (a[2] * a[1]) / square / 4.0;

        K[2][0] = (a[0] * a[2]) / square / 4.0;
        K[2][1] = (a[1] * a[2]) / square / 4.0;
        K[2][2] = (a[2] * a[2]) / square / 4.0;
    }

    void mass(double ** K, double square)
    {
        K[0][0] = K[1][1] = K[2][2] = square / 6.;
        K[0][1] = K[0][2] = K[1][0] = K[1][2] = K[2][0] = K[2][1] = square / 12.;
    }

    void testP(double square, double *y, double *rho, double *vector)
    {
        vector[0] = 9.81 * square * (rho[0] * (1.0 - y[0]) / 6.0 + rho[1] * (1.0 - y[1]) / 12.0 + rho[2] *(1.0 - y[2]) / 12.0);
        vector[1] = 9.81 * square * (rho[0] * (1.0 - y[0]) / 12.0 + rho[1] *(1.0 - y[1]) / 6.0 + rho[2] * (1.0 - y[2]) / 12.0);
        vector[2] = 9.81 * square * (rho[0] * (1.0 - y[0]) / 12.0 + rho[1] *(1.0 - y[1]) / 12.0 + rho[2] *(1.0 - y[2]) / 6.0);
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
};

#endif //NEWFEM_LOCALMATRIX_H

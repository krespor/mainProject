//
// Created by konstantin on 28.02.19.
//

#include "SUPG.h"

SUPG::SUPG(Mesh mesh, Arguments arguments, string pathToResult) : Solver(mesh, arguments, pathToResult)
{
    reservMemory();
    init();

    calc();

}

void SUPG::reservMemory()
{
    c = new double[mesh.n];
    cn = new double[mesh.n];
    u = new double[mesh.n];
    v = new double[mesh.n];

    localU = new double[3];
    localV = new double[3];

    localVector0 = new double[3];
    localVector1 = new double[3];
    localVector2 = new double[3];
    localVector3 = new double[3];

    localMatrix0 = new double*[3];
    localMatrix1 = new double*[3];
    localMatrix2 = new double*[3];
    localMatrix3 = new double*[3];

    for (int i = 0; i < 3; i++)
    {
        localMatrix0[i] = new double[3];
        localMatrix1[i] = new double[3];
        localMatrix2[i] = new double[3];
        localMatrix3[i] = new double[3];
    }
}

void SUPG::init()
{
    methods.null(c, mesh.n);
    methods.null(cn, mesh.n);
    //methods.null(v, mesh.n);

    methods.actionsVC(u, 1.0 / sqrt(2), mesh.n, '=');
    methods.actionsVC(v, 1.0 / sqrt(2), mesh.n, '=');

    hElem = 0;
}

void SUPG::fillSLAE()
{
    for (int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];

            localVector0[i] = cn[k[i]];
            localU[i] = u[k[i]];
            localV[i] = v[k[i]];

            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        calc_a_b();

        hElem = calcH(localU, localV, mesh.square[t]);

        supgMatrixMass(a, b, localU, localV, hElem, mesh.square[t], 0, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        //supgMatrix(b, mesh.square[t], hElem, localMatrix1);
        supgFull(a, b, localU, localV, hElem, mesh.square[t], 0, localMatrix1);

        supgMatrixMass(a, b, localU, localV, hElem, mesh.square[t], 0, localMatrix2);
        methods.multMV(localMatrix2, localVector0, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j];
            }
        }
    }
}

void SUPG::calc()
{
    /*for (int i = 0; i < mesh.n; i++)
        if ((mesh.nodes[i][0] >= 1.85) && (mesh.nodes[i][0] <= 2.15))
            if ((mesh.nodes[i][1] >= 0.35) && (mesh.nodes[i][1] <= 0.65))
                cn[i] = 1;
*/

    cn[0] = 1;
    cn[1] = 1;
    cn[2] = 1;
    cn[3] = 1;
    cn[225] = 1;
    cn[224] = 1;
    cn[226] = 1;
    cn[227] = 1;
    cn[228] = 1;
    cn[229] = 1;
    cn[230] = 1;
    cn[231] = 1;
    cn[1310] = 1;
    cn[1311] = 1;
    cn[1308] = 1;
    cn[1309] = 1;
    cn[1307] = 1;

    recordData(vector<double *>{cn}, vector<string>{"c"});

    while (runTime < endTime)
    {
        runTime = del_t * countIterations;
        cout << "time = " << runTime << endl;
        fillSLAE();

        solveSLAE(c);

        autoRecordData(vector<double *>{c}, vector<string>{"c"});

        methods.equateV(cn, c, mesh.n);

        countIterations++;
    }
}

double SUPG::calcH(double *u, double *v, double square)
{
    /*double inf = std::numeric_limits<double>::infinity();
    double *lineOne = new double [4];
    double *lineTwo = new double [4];
    int temp = 0;
    double lengthLineOne;

    double x1, x2, x3, x4, y1, y2, y3, y4;
    double k1, k2, b1, b2;
    double hx, hy;

    double v1 = (u[0] + u[1] + u[2]) / 3.0;
    double v2 = (v[0] + v[1] + v[2]) / 3.0;

    if ((fabs(v1) < 0.00000001) && (fabs(v2) < 0.00000001))
    {
        return 0.0;
    }

    for (int i = 0; i<3; i++)
    {
        lineOne[0] = x[i];
        lineOne[1] = y[i];
        lineOne[2] = x[i] + v1 * 10;
        lineOne[3] = y[i] + v2 * 10;

        lengthLineOne = methods.findLength(lineOne[0], lineOne[2], lineOne[1], lineOne[3]);

        for (int j = 900; lengthLineOne < 2.0; j *=900)
        {
            lineOne[2] = x[i] + v1 * j;
            lineOne[3] = y[i] + v2 * j;

            lengthLineOne = methods.findLength(lineOne[0], lineOne[2], lineOne[1], lineOne[3]);
        }

        for(int j = 0; j < 3; j++)
        {
            if (i != j)
            {
                lineTwo[temp] = x[j];
                lineTwo[temp + 1] = y[j];

                temp += 2;
            }
        }

        temp = 0;

        x1 = lineOne[0];
        y1 = lineOne[1];
        x2 = lineOne[2];
        y2 = lineOne[3];

        x3 = lineTwo[0];
        y3 = lineTwo[1];
        x4 = lineTwo[2];
        y4 = lineTwo[3];

        double tmp;

        if (x1 >= x2)
        {
            tmp = x1;
            x1 = x2;
            x2 = tmp;

            tmp = y1;
            y1 = y2;
            y2 = tmp;
        }

        if (x3 >= x4)
        {
            tmp = x3;
            x3 = x4;
            x4 = tmp;

            tmp = y3;
            y3 = y4;
            y4 = tmp;
        }

        k1 =  ( y2 - y1 ) / ( x2 - x1 );
        k2 =  ( y4 - y3 ) / ( x4 - x3 );

        if (k1 == inf)
        {
            k1 = 999999999;
        }
        if (k2 == inf)
        {
            k2 = 999999999;
        }
        if (k1 == -inf)
        {
            k1 = -999999999;
        }
        if (k2 == -inf)
        {
            k2 = -999999999;
        }

        b1 = y1 - k1 * x1;
        b2 = y3  - k2 * x3;

        hx = ( b2 - b1 ) / ( k1 -  k2 );
        hy = k1 * hx + b1;

        tmp = fabs(methods.findLength(lineTwo[0],lineTwo[2],lineTwo[1],lineTwo[3]) - (methods.findLength(lineTwo[0],hx,lineTwo[1], hy) + methods.findLength(hx,lineTwo[2] ,hy,lineTwo[3])));
        if (tmp < 10e-6)
        {
            tmp = fabs(methods.findLength(lineOne[0],lineOne[2],lineOne[1],lineOne[3]) - (methods.findLength(lineOne[0],hx,lineOne[1], hy) + methods.findLength(hx, lineOne[2],hy,lineOne[3])));
            if (tmp < 10e-6)
            {
                delete lineOne;
                delete lineTwo;
                return methods.findLength(x[i], hx, y[i], hy);
            }
        }

    }

    delete lineOne;
    delete lineTwo;

    return  0.0;*/



    double uAvg = (u[0] + u[1] + u[2]) / 3.0;
    double vAvg = (v[0] + v[1] + v[2]) / 3.0;

    double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

    double he = 2.0 / (
            fabs((uAvg * b[0]) / (2.0 * square * uMod) + (vAvg * a[0]) / (2.0 * square * uMod)) +
            fabs((uAvg * b[1]) / (2.0 * square * uMod) + (vAvg * a[1]) / (2.0 * square * uMod)) +
            fabs((uAvg * b[2]) / (2.0 * square * uMod) + (vAvg * a[2]) / (2.0 * square * uMod)));

    return he;
}

void SUPG::supgMatrixMass(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix)
{
        double uAvg = (u[0] + u[1] + u[2]) / 3.0;
        double vAvg = (v[0] + v[1] + v[2]) / 3.0;
        double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

        if (uMod != 0)
        {
            double alpha = 1;//(1.0 / tanh((uMod * h) / (2 * k))) - (2.0 * k) / (uMod * h);

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

void SUPG::supgMatrix(double *b, double square, double h, double **matrix)
{
    matrix[0][0] = b[0] / 6.0 /**/ + (b[0] * b[0] * h) / (8.0 * square);
    matrix[0][1] = b[1] / 6.0 /**/ + (b[1] * b[0] * h) / (8.0 * square);
    matrix[0][2] = b[2] / 6.0 /**/ + (b[2] * b[0] * h) / (8.0 * square);

    matrix[1][0] = b[0] / 6.0 /**/ + (b[0] * b[1] * h) / (8.0 * square);
    matrix[1][1] = b[1] / 6.0 /**/ + (b[1] * b[1] * h) / (8.0 * square);
    matrix[1][2] = b[2] / 6.0 /**/ + (b[2] * b[1] * h) / (8.0 * square);

    matrix[2][0] = b[0] / 6.0 /**/ + (b[0] * b[2] * h) / (8.0 * square);
    matrix[2][1] = b[1] / 6.0 /**/ + (b[1] * b[2] * h) / (8.0 * square);
    matrix[2][2] = b[2] / 6.0 /**/ + (b[2] * b[2] * h) / (8.0 * square);
}

void SUPG::supgFull(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix)
{
    if (true)
    {
        double uAvg = (u[0] + u[1] + u[2]) / 3.0;
        double vAvg = (v[0] + v[1] + v[2]) / 3.0;
        double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

        if (uMod != 0)
        {
            double alpha = 1;//(1.0 / tanh((uMod * h) / (2 * k))) - (2.0 * k) / (uMod * h);

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
            methods.null(matrix, 3);
        }

    }

}

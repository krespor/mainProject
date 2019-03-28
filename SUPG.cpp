//
// Created by konstantin on 28.02.19.
//

#include "SUPG.h"

SUPG::SUPG(Mesh mesh, Arguments arguments, string pathToResult) : Solver(mesh, arguments, pathToResult)
{
    d = arguments.individualArguments["D"];
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

    methods.actionsVC(u, 1, mesh.n, '=');
    methods.actionsVC(v, 0, mesh.n, '=');

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

        //localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], d, localMatrix0);
        localMatrix.mass(localMatrix0, mesh.square[t]);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.convectiveMembers(localMatrix1, localU, localV, a, b);
        //localMatrix.supg.convectiveMembers(a, b, localU, localV, hElem, mesh.square[t], d, localMatrix1);

        localMatrix.mass(localMatrix2, mesh.square[t]);
        //localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], d, localMatrix2);
        methods.multMV(localMatrix2, localVector0, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

        supgMatrixLaplas(localU, localV, hElem, mesh.square[t], d, localMatrix2);
        methods.multMC(localMatrix2, d, 3);

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix2[i][j];
            }
        }
    }
}

void SUPG::calc()
{
    /*for (int i = 0; i < mesh.n; i++)
        if ((mesh.nodes[i][0] >= 1.85) && (mesh.nodes[i][0] <= 2.15))
            if ((mesh.nodes[i][1] >= 0.35) && (mesh.nodes[i][1] <= 0.65))
                cn[i] = 1;*/

    /*cn[0] = 1;
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
    cn[1307] = 1;*/

    int tb;
    for (int i = 0; i < mesh.border[0].n;i++)
    {
        tb = mesh.border[0].numberNodes[i];

        if (mesh.nodes[tb][1] <= 0.65 && mesh.nodes[tb][1] >= 0.35) {
            cn[tb] = 1;
        }

    }

    recordData(vector<double *>{cn}, vector<string>{"c"});

    while (runTime < endTime)
    {
        runTime = del_t * countIterations;
        cout << "time = " << runTime << endl;
        fillSLAE();
        //conditionBorder_1(0, 0);
        //conditionBorder_1(0, 1);
        //conditionBorder_1(0, 2);
        //conditionBorder_1(0, 3);
        solveSLAE(c);

        for (int i = 0; i < mesh.border[0].n;i++)
        {
            tb = mesh.border[0].numberNodes[i];

            if (mesh.nodes[tb][1] <= 0.65 && mesh.nodes[tb][1] >= 0.35) {
                c[tb] = 1;
            }

        }

        /*for (int i = 0; i < mesh.n; i++)
            if ((mesh.nodes[i][0] >= 1.85) && (mesh.nodes[i][0] <= 2.15))
                if ((mesh.nodes[i][1] >= 0.35) && (mesh.nodes[i][1] <= 0.65))
                    c[i] = 1;*/

        autoRecordData(vector<double *>{c}, vector<string>{"c"});

        methods.equateV(cn, c, mesh.n);

        countIterations++;
    }
}

double SUPG::calcH(double *u, double *v, double square)
{
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
            double alpha = (1.0 / tanh((uMod * h) / (2 * k))) - (2.0 * k) / (uMod * h);

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

void SUPG::supgFull(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix)
{
    if (true)
    {
        double uAvg = (u[0] + u[1] + u[2]) / 3.0;
        double vAvg = (v[0] + v[1] + v[2]) / 3.0;
        double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

        if (uMod != 0)
        {
            double alpha = (1.0 / tanh((uMod * h) / (2 * k))) - (2.0 * k) / (uMod * h);

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

void SUPG::supgMatrixLaplas(double *u, double *v, double h, double square, double k, double **matrix)
{
    double uAvg = (u[0] + u[1] + u[2]) / 3.0;
    double vAvg = (v[0] + v[1] + v[2]) / 3.0;
    double uMod = sqrt(uAvg * uAvg + vAvg * vAvg);

    double pe = (uMod * h) / (2.0 * k);
    double tau = (h / (2.0 * uMod)) *
            (1.0 / tanh(pe) - 1.0 / pe);

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
}

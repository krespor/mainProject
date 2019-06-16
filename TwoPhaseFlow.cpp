//
// Created by konstantin on 06.04.19.
//

#include "TwoPhaseFlow.h"

TwoPhaseFlow::TwoPhaseFlow(Mesh mesh, Arguments arguments, string pathToResult) : Solver(mesh, arguments,pathToResult)
{
    rho0 = arguments.individualArguments["rho0"];
    rho1 = arguments.individualArguments["rho1"];
    mu0 = arguments.individualArguments["mu0"];
    mu1 = arguments.individualArguments["mu1"];

    reservMemory();
    init();

    calc();
}

void TwoPhaseFlow::reservMemory()
{
    p = new double[mesh.n];
    u = new double[mesh.n];
    v = new double[mesh.n];
    un = new double[mesh.n];
    vn = new double[mesh.n];
    uStar = new double[mesh.n];
    vStar = new double[mesh.n];

    alpha = new double[mesh.n];
    alpha_n = new double[mesh.n];

    mu = new double[mesh.n];
    rho = new double[mesh.n];

    localU = new double[3];
    localV = new double[3];
    localP = new double[3];
    localRho = new double[3];
    localMu = new double[3];
    localAlpha_n = new double[3];

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

void TwoPhaseFlow::init()
{
    g = 9.81;
    methods.null(p, mesh.n);
    methods.null(uStar, mesh.n);
    methods.null(vStar, mesh.n);
    methods.null(u, mesh.n);
    methods.null(v, mesh.n);
    methods.null(un, mesh.n);
    methods.null(vn, mesh.n);

    methods.null(alpha, mesh.n);
    methods.null(alpha_n, mesh.n);

    methods.null(mu, mesh.n);
    methods.null(rho, mesh.n);
}

void TwoPhaseFlow::calc()
{
    for (int i = 0; i < mesh.n; i++)
    {
        if (mesh.nodes[i][1] >= 0.5)
            alpha_n[i] = alpha[i] = 1;
        //if ((mesh.nodes[i][1] <= 0.8) && (mesh.nodes[i][0] <= 0.3))
        //    alpha_n[i] = alpha[i] = 1;
        //if ((mesh.nodes[i][1] <= 0.7) && (mesh.nodes[i][0] >= 0.4) && (mesh.nodes[i][0] <= 0.6))
        //    alpha_n[i] = alpha[i] = 1;

        //if ((mesh.nodes[i][1] <= 0.18) && (mesh.nodes[i][1] >= 0.12) && (mesh.nodes[i][0] >= 0.47) && (mesh.nodes[i][0] <= 0.53))
        //    alpha_n[i] = alpha[i] = 0;
    }

    setPhase();

    recordData(vector<double *>{uStar, vStar, p, u, v, mu, rho, alpha}, vector<string>{"u*", "v*", "p", "u", "v", "mu", "rho", "alpha"});

    while (runTime < endTime)
    {
        runTime = del_t * countIterations;

        cout << "U*" << endl;
        fillSLAE_uStar();
        solveSLAE(uStar);

        cout << "V*" << endl;
        fillSLAE_vStar();
        solveSLAE(vStar);

        cout << "P" << endl;
        fillSLAE_p();
        bcPressure();
        //conditionBorder_2(g, 1);
        conditionBorder_1(0, 0);
        solveSLAE(p);

        cout << "U" << endl;
        fillSLAE_u();
        solveSLAE(u);

        cout << "V" << endl;
        fillSLAE_v();
        solveSLAE(v);

        currectUV();

        cout << "alpha" << endl;
        fillSLAE_alpha();
        //conditionBorder_1(0, 0);
        solveSLAE(alpha);

        for (int i = 0; i < mesh.n; i++)
        {
            if (alpha[i] < 0)
                alpha[i] = 0;
            else if (alpha[i] > 1)
                alpha[i] = 1;
        }

        setPhase();

        methods.equateV(alpha_n, alpha, mesh.n);
        methods.equateV(un, u, mesh.n);
        methods.equateV(vn, v, mesh.n);

        autoRecordData(vector<double *>{uStar, vStar, p, u, v, mu, rho, alpha}, vector<string>{"u*", "v*", "p", "u", "v", "mu", "rho", "alpha"});
        countIterations++;
    }
}

void TwoPhaseFlow::fillSLAE_uStar()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = un[k[i]];
            localV[i] = vn[k[i]];
            localRho[i] = rho[k[i]];
            localMu[i] = mu[k[i]];
        }
        calc_a_b();

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.twoPhase.convectiveMembersRho(a, b, localU, localV, localRho, mesh.square[t], localMatrix1);

        localMatrix.laplass(localMatrix2, mesh.square[t], a, b);
        methods.multMC(localMatrix2, mu1, 3);
        //localMatrix.twoPhase.laplassMuX(a, b, localMu, localV, mesh.square[t], localVector1, localMatrix2);

        localMatrix.mass(localMatrix3, mesh.square[t]);
        methods.multMV(localMatrix3, localU, localVector0, 3);
        methods.actionsVC(localVector0, 1.0 / del_t, 3, '*');


        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i];// + localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix2[i][j];
            }
        }
    }
}

void TwoPhaseFlow::fillSLAE_vStar()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = un[k[i]];
            localV[i] = vn[k[i]];
            localRho[i] = rho[k[i]];
            localMu[i] = mu[k[i]];
        }
        calc_a_b();

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.twoPhase.convectiveMembersRho(a, b, localU, localV, localRho, mesh.square[t], localMatrix1);

        localMatrix.laplass(localMatrix2, mesh.square[t], a, b);
        methods.multMC(localMatrix2, mu1, 3);
        //localMatrix.twoPhase.laplassMuY(a, b, localMu, localU, mesh.square[t], localVector2, localMatrix2);

        localMatrix.mass(localMatrix3, mesh.square[t]);
        methods.multMV(localMatrix3, localRho, localVector0, 3);
        methods.actionsVC(localVector0, -g, 3, '*');

        localMatrix.mass(localMatrix3, mesh.square[t]);
        methods.multMV(localMatrix3, localV, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i] + localVector1[i];// + localVector2[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix2[i][j];
            }
        }
    }
}

void TwoPhaseFlow::fillSLAE_p()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = uStar[k[i]];
            localV[i] = vStar[k[i]];
            localRho[i] = rho[k[i]];
        }
        calc_a_b();

        localMatrix.twoPhase.laplassP(a, b, localRho, mesh.square[t], localMatrix0);

        localMatrix.derivative(localMatrix1, a);
        methods.multMV(localMatrix1, localV, localVector0, 3);
        methods.actionsVC(localVector0, 1.0 / del_t, 3, '*');

        localMatrix.derivative(localMatrix1, a);
        methods.multMV(localMatrix1, localU, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i] + localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += -localMatrix0[i][j];
            }
        }
    }
}

void TwoPhaseFlow::fillSLAE_u()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = uStar[k[i]];
            localP[i] = p[k[i]];
            localRho[i] = rho[k[i]];
        }
        calc_a_b();

        localMatrix.mass(localMatrix0, mesh.square[t]);

        localMatrix.twoPhase.derPPP(b, localRho, localP, localVector0);
        methods.actionsVC(localVector0, -del_t, 3, '*');

        localMatrix.mass(localMatrix1, mesh.square[t]);
        methods.multMV(localMatrix1, localU, localVector1, 3);

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i] + localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j];
            }
        }
    }
}

void TwoPhaseFlow::fillSLAE_v()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localV[i] = vStar[k[i]];
            localP[i] = p[k[i]];
            localRho[i] = rho[k[i]];
        }
        calc_a_b();

        localMatrix.mass(localMatrix0, mesh.square[t]);

        localMatrix.twoPhase.derPPP(a, localRho, localP, localVector0);
        methods.actionsVC(localVector0, -del_t, 3, '*');

        localMatrix.mass(localMatrix1, mesh.square[t]);
        methods.multMV(localMatrix1, localV, localVector1, 3);

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i] + localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j];
            }
        }
    }
}

void TwoPhaseFlow::fillSLAE_alpha()
{
    double hElem;
    for (int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];

            localAlpha_n[i] = alpha_n[k[i]];
            localU[i] = u[k[i]];
            localV[i] = v[k[i]];
        }

        calc_a_b();

        hElem = calcH(localU, localV, mesh.square[t]);

        localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], 0, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.supg.convectiveMembers(a, b, localU, localV, hElem, mesh.square[t], 0, localMatrix1);

        localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], 0, localMatrix2);
        methods.multMV(localMatrix2, localAlpha_n, localVector1, 3);
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

void TwoPhaseFlow::currectUV()
{
    int t;
    for (int j = 0; j < mesh.countBorder; j++)
        for (int i = 0; i<mesh.border[j].n; i++)
        {
            t = mesh.border[j].numberNodes[i];
            u[t] = v[t] = 0;
        }
}

double TwoPhaseFlow::calcH(double *u, double *v, double square)
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

void TwoPhaseFlow::setPhase()
{
    for (int i = 0; i<mesh.n; i++)
    {
        rho[i] = (1 - alpha[i]) * rho0 + alpha[i] * rho1;
        mu[i] = (1 - alpha[i]) * mu0 + alpha[i] * mu1;
    }
}

void TwoPhaseFlow::bcPressure()
{
    unsigned int numberBorder = 1;

    double f;  // среднее значение функции на отрезке
    double X1, Y1, X2, Y2;   //координаты отрезка
    double l; //длинна отрезка
    int k1, k2;   //номера точек отрезка

    for (unsigned int k = 0; k < mesh.border[numberBorder].m; k++) //цикл по отрезку// {
    {
        k1 = mesh.border[numberBorder].line[k][0];     //
        k2 = mesh.border[numberBorder].line[k][1]; //  находим номера точек первого отрезка

        X1 = mesh.nodes[k1][0]; //
        Y1 = mesh.nodes[k1][1]; // берём координаты первого отрезка
        X2 = mesh.nodes[k2][0]; //
        Y2 = mesh.nodes[k2][1]; //

        double c1 = -vStar[k1] / del_t;
        double c2 = -vStar[k2] / del_t;

        l = methods.findLength(X1, X2, Y1, Y2); //нахождение длинны отрезка
        f = (c1 + c2) / 2; //нахождение среднего значения функции на отрезке
        B[k1] -= l / 2 * f;
        B[k2] -= l / 2 * f;
    }

    numberBorder = 2;

    for (unsigned int k = 0; k < mesh.border[numberBorder].m; k++) //цикл по отрезку// {
    {
        k1 = mesh.border[numberBorder].line[k][0];     //
        k2 = mesh.border[numberBorder].line[k][1]; //  находим номера точек первого отрезка

        X1 = mesh.nodes[k1][0]; //
        Y1 = mesh.nodes[k1][1]; // берём координаты первого отрезка
        X2 = mesh.nodes[k2][0]; //
        Y2 = mesh.nodes[k2][1]; //

        double c1 = uStar[k1] / del_t;
        double c2 = uStar[k2] / del_t;

        if (X1 == 0)
        {
            c1 *= -1;
            c2 *= -1;
        }

        l = methods.findLength(X1, X2, Y1, Y2); //нахождение длинны отрезка
        f = (c1 + c2) / 2; //нахождение среднего значения функции на отрезке
        B[k1] -= l / 2 * f;
        B[k2] -= l / 2 * f;
    }
}








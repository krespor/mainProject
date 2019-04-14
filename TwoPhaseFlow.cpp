//
// Created by konstantin on 06.04.19.
//

#include "TwoPhaseFlow.h"

TwoPhaseFlow::TwoPhaseFlow(Mesh mesh, Arguments arguments, string pathToResult) : Solver(mesh, arguments,pathToResult)
{
    rho1 = arguments.individualArguments["rho1"];
    rho2 = arguments.individualArguments["rho2"];
    mu1 = arguments.individualArguments["mu1"];
    mu2 = arguments.individualArguments["mu2"];

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
        if (mesh.nodes[i][1] <= 0.5)
        {
            alpha_n[i] = alpha[1] = 1;
        }
    }

    setPhase();

    recordData(vector<double *>{u, v, p, alpha_n, mu, rho}, vector<string>{"u", "v", "p", "alpha", "mu", "rho"});

    while (runTime < endTime)
    {
        runTime = del_t * countIterations;
        cout << "time = " << runTime << endl;
        setPhase();

        /*fillSLAE_uStar();
        conditionBorder_1(0, 1);
        solveSLAE(uStar);

        fillSLAE_vStar();
        conditionBorder_1(0, 1);
        solveSLAE(vStar);*/

        fillSLAE_p();
        conditionBorder_1(0, 0);
        solveSLAE(p);

        /*fillSLAE_u();
        solveSLAE(u);

        fillSLAE_v();
        solveSLAE(v);

        currectUV();*/

        /*fillSLAE_alpha();
        conditionBorder_1(0, 0);
        solveSLAE(alpha);

        methods.equateV(alpha_n, alpha, mesh.n);*/
        methods.equateV(un, u, mesh.n);
        methods.equateV(vn, v, mesh.n);

        autoRecordData(vector<double *>{uStar, vStar, p, alpha, mu, rho}, vector<string>{"u", "v", "p", "alpha", "mu", "rho"});
        countIterations++;
        cout << "A" << endl;
        cin.get();
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

        localMatrix.twoPhase.laplassMu(a, b, localMu, mesh.square[t], localMatrix2);

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix3);
        methods.multMV(localMatrix3, localU, localVector0, 3);
        methods.actionsVC(localVector0, 1.0 / del_t, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix2[i][j];  //
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

        localMatrix.twoPhase.laplassMu(a, b, localMu, mesh.square[t], localMatrix2);

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix3);
        methods.multMV(localMatrix3, localV, localVector0, 3);
        methods.actionsVC(localVector0, 1.0 / del_t, 3, '*');

        localMatrix.mass(localMatrix3, mesh.square[t]);
        methods.multMV(localMatrix3, localRho, localVector1, 3);
        methods.actionsVC(localVector1, -g, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i] + localVector1[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix2[i][j];  //
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

        localMatrix.laplass(localMatrix0, mesh.square[t], a, b);

        //localMatrix.twoPhase.derivativeRho(a, b, localRho, localU, localV, localVector0);
        //methods.actionsVC(localVector0, 1.0/del_t, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j];
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

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.derivative(localMatrix1, b);
        methods.multMV(localMatrix1, localP, localVector0, 3);
        methods.actionsVC(localVector0, -1, 3, '*');

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix1);
        methods.multMV(localMatrix1, localU, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

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

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.derivative(localMatrix1, b);
        methods.multMV(localMatrix1, localP, localVector0, 3);
        methods.actionsVC(localVector0, -1, 3, '*');

        localMatrix.twoPhase.massRho(mesh.square[t], localRho, localMatrix1);
        methods.multMV(localMatrix1, localV, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

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

    for (int i = 0; i<mesh.border[1].n; i++)
    {
        t = mesh.border[1].numberNodes[i];
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
        if (alpha_n[i] > 0.98)
        {
            mu[i] = mu2;
            rho[i] = rho2;
        } else
        {
            mu[i] = mu1;
            rho[i] = rho1;
        }
    }
}

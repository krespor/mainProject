//
// Created by konstantin on 06.03.19.
//

#include "NavierStokes.h"

NavierStokes::NavierStokes() : Solver()
{
    rho = arguments.individualArguments["rho"];
    mu = arguments.individualArguments["mu"];

    reservMemory();
    init();

    calcBackwardStep();
}

void NavierStokes::reservMemory()
{
    p = new double[mesh.n];
    u = new double[mesh.n];
    v = new double[mesh.n];
    un = new double[mesh.n];
    vn = new double[mesh.n];
    uStar = new double[mesh.n];
    vStar = new double[mesh.n];

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

void NavierStokes::init()
{
    methods.null(p, mesh.n);
    methods.null(uStar, mesh.n);
    methods.null(vStar, mesh.n);
    methods.null(u, mesh.n);
    methods.null(v, mesh.n);
    methods.null(un, mesh.n);
    methods.null(vn, mesh.n);
}

void NavierStokes::calcBackwardStep()
{
    recordData(vector<double *>{u, v}, vector<string>{"u", "v"});

    while (runTime < endTime)
    {
        runTime = del_t * countIterations;

        fillSLAE_uStar();
        privateBorderCondition();
        conditionBorder_1(0, 1);
        solveSLAE(uStar);

        fillSLAE_vStar();
        conditionBorder_1(0, 0);
        conditionBorder_1(0, 1);
        solveSLAE(vStar);

        fillSLAE_p();
        conditionBorder_2(12 * mu, 0);
        conditionBorder_1(0, 2);
        solveSLAE(p);

        fillSLAE_u();
        solveSLAE(u);

        fillSLAE_v();
        solveSLAE(v);

        currectUV();

        methods.equateV(un, u, mesh.n);
        methods.equateV(vn, v, mesh.n);

        autoRecordData(vector<double *>{u, v}, vector<string>{"u", "v"});
        countIterations++;
    }
}

void NavierStokes::fillSLAE_uStar()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = un[k[i]];
            localV[i] = vn[k[i]];
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        calc_a_b();

        localMatrix.mass(localMatrix0, mesh.square[t]);
        methods.multMC(localMatrix0, 1. / del_t, 3);

        localMatrix.laplass(localMatrix1, mesh.square[t], a, b);
        methods.multMC(localMatrix1, mu / rho, 3);

        localMatrix.mass(localMatrix2, mesh.square[t]);
        methods.multMV(localMatrix2, localU, localVector0, 3);
        methods.actionsVC(localVector0, 1. / del_t, 3, '*');

        localMatrix.convectiveMembers(localMatrix3, localU, localV, a, b);

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix3[i][j];
            }
        }
    }
}

void NavierStokes::fillSLAE_vStar()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = un[k[i]];
            localV[i] = vn[k[i]];
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        calc_a_b();

        localMatrix.mass(localMatrix0, mesh.square[t]);
        methods.multMC(localMatrix0, 1. / del_t, 3);

        localMatrix.laplass(localMatrix1, mesh.square[t], a, b);
        methods.multMC(localMatrix1, mu / rho, 3);

        localMatrix.mass(localMatrix2, mesh.square[t]);
        methods.multMV(localMatrix2, localV, localVector0, 3);
        methods.actionsVC(localVector0, 1. / del_t, 3, '*');

        localMatrix.convectiveMembers(localMatrix3, localU, localV, a, b);

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i];
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += localMatrix0[i][j] + localMatrix1[i][j] + localMatrix3[i][j];  //
            }
        }
    }
}

void NavierStokes::fillSLAE_p()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = uStar[k[i]];
            localV[i] = vStar[k[i]];
        }
        for (unsigned int i = 0; i<3; i++)
        {
            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        calc_a_b();

        localMatrix.laplass(localMatrix0, mesh.square[t], a, b);

        localMatrix.derivative(localMatrix2, b);
        methods.multMV(localMatrix2, localU, localVector0, 3);
        methods.actionsVC(localVector0, 1.0 / del_t, 3, '*');

        localMatrix.derivative(localMatrix2, a);
        methods.multMV(localMatrix2, localV, localVector1, 3);
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

void NavierStokes::fillSLAE_u()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = p[k[i]];
            localV[i] = uStar[k[i]];
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        calc_a_b();

        localMatrix.mass(localMatrix0, mesh.square[t]);

        localMatrix.derivative(localMatrix1, b);
        methods.multMV(localMatrix1, localU, localVector0, 3);
        methods.actionsVC(localVector0, -del_t / rho, 3, '*');

        localMatrix.mass(localMatrix2, mesh.square[t]);
        methods.multMV(localMatrix2, localV, localVector1, 3);


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

void NavierStokes::fillSLAE_v()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
            localU[i] = p[k[i]];
            localV[i] = vStar[k[i]];
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        calc_a_b();

        localMatrix.mass(localMatrix0, mesh.square[t]);

        localMatrix.derivative(localMatrix1, a);
        methods.multMV(localMatrix1, localU, localVector0, 3);
        methods.actionsVC(localVector0, -del_t / rho, 3, '*');

        localMatrix.mass(localMatrix2, mesh.square[t]);
        methods.multMV(localMatrix2, localV, localVector1, 3);

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

void NavierStokes::currectUV()
{
    double y;
    int t;
    int nubmerBorder = 1;

    for (int i = 0; i < mesh.border[nubmerBorder].n; i++)
    {
        t = mesh.border[nubmerBorder].numberNodes[i];
        u[t] = v[t] = 0;
    }

    nubmerBorder = 0;
    for (int i = 0; i < mesh.border[nubmerBorder].n; i++)
    {
        t = mesh.border[nubmerBorder].numberNodes[i];
        y = mesh.nodes[t][1];
        v[t] = 0;
        u[t] = 6 * (y - y * y);
    }
}

void NavierStokes::privateBorderCondition()
{
    unsigned int numberBorder = 0; //номер границы
    unsigned int t; //точка на границе
    unsigned int countDotBorder; //количество точек на границе
    double y, c;

    countDotBorder = mesh.border[numberBorder].n;

    for (unsigned int i = 0; i < countDotBorder; i++)  //цикл по границе
    {
        t = mesh.border[numberBorder].numberNodes[i]; //выбираем точку на границе
        y = mesh.nodes[t][1];
        c = 6 * (y - y * y);
        conditionNode_1(c, t);
    }
}
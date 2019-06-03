//
// Created by konstantin on 06.03.19.
//

#include "NavierStokes.h"

NavierStokes::NavierStokes(Mesh mesh, Arguments arguments, string pathToResult) : Solver(mesh, arguments,pathToResult)
{
    rho = arguments.individualArguments["rho"];
    mu = arguments.individualArguments["mu"];
    d = arguments.individualArguments["D"];

    reservMemory();
    init();

    //calcCaverna();
    calcTwoPhaseFlow();
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

    c = new double[mesh.n];
    cn = new double[mesh.n];

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

    methods.null(c, mesh.n);
    methods.null(cn, mesh.n);
}

void NavierStokes::calcBackwardStep()
{
    int tb;

    recordData(vector<double *>{u, v}, vector<string>{"u", "v"});

    while (runTime < endTime)
    {
        runTime = del_t * countIterations;
        cout << "time = " << runTime << endl;

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

        /*fillSLAE_c();
        for (int i = 0; i < mesh.border[0].n;i++)
        {
            tb = mesh.border[0].numberNodes[i];

            if (mesh.nodes[tb][1] <= 0.65 && mesh.nodes[tb][1] >= 0.35) {
                conditionNode_1(1, tb);
            } else
            {
                conditionNode_1(0, tb);
            }
        }
        solveSLAE(c);

        methods.equateV(cn, c, mesh.n);*/
        methods.equateV(un, u, mesh.n);
        methods.equateV(vn, v, mesh.n);

        autoRecordData(vector<double *>{u, v}, vector<string>{"u", "v"});
        countIterations++;
    }
}


void NavierStokes::calcTwoPhaseFlow()
{
    while (runTime < endTime) {
        runTime = del_t * countIterations;
        cout << "time = " << runTime << endl;

        fillSLAE_uStar();
        solveSLAE(uStar);

        fillSLAE_vStar();
        solveSLAE(vStar);

        fillSLAE_p();
        conditionBorder_2(1000*9.81, 1);
        conditionBorder_1(0, 0);

        solveSLAE(p);

        fillSLAE_u();
        solveSLAE(u);

        fillSLAE_v();
        solveSLAE(v);

        for (int i = 0; i < mesh.countBorder; i++)
            for (int j = 0; j < mesh.border[i].n; j++)
            {
                int tb = mesh.border[i].numberNodes[j];
                v[tb] = u[tb] = 0;
            }

        methods.equateV(un, u, mesh.n);
        methods.equateV(vn, v, mesh.n);

        recordData(vector<double *>{uStar, vStar, p, u, v}, vector<string>{"u*", "v*", "p", "u", "v"});
        countIterations++;
    }
}



void NavierStokes::fillSLAE_uStar()
{
    //double hElem;
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
        //hElem = calcH(localU, localV, mesh.square[t]);

        localMatrix.mass(localMatrix0, mesh.square[t]);
        //localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix0);
        methods.multMC(localMatrix0, 1. / del_t, 3);

        localMatrix.laplass(localMatrix1, mesh.square[t], a, b);
        //localMatrix.supg.laplass(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix1);
        methods.multMC(localMatrix1, mu / rho, 3);

        localMatrix.mass(localMatrix2, mesh.square[t]);
        //localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix2);
        methods.multMV(localMatrix2, localU, localVector0, 3);
        methods.actionsVC(localVector0, 1. / del_t, 3, '*');

        /*for (int i = 0; i < 3; i++)
        {
            if (x[i] == 53)
            {
                localU[i] = 0;
                localV[i] = 0;
            }
        }*/

        localMatrix.convectiveMembers(localMatrix3, localU, localV, a, b);
        //localMatrix.supg.convectiveMembers(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix3);

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

void NavierStokes::fillSLAE_vStar()
{
    //double hElem;
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
        //hElem = calcH(localU, localV, mesh.square[t]);

        localMatrix.mass(localMatrix0, mesh.square[t]);
        //localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix0);
        methods.multMC(localMatrix0, 1. / del_t, 3);

        localMatrix.laplass(localMatrix1, mesh.square[t], a, b);
        //localMatrix.supg.laplass(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix1);
        methods.multMC(localMatrix1, mu / rho, 3);

        localMatrix.mass(localMatrix2, mesh.square[t]);
        //localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix2);
        methods.multMV(localMatrix2, localV, localVector0, 3);
        methods.actionsVC(localVector0, 1. / del_t, 3, '*');

        /*for (int i = 0; i < 3; i++)
        {
            if (x[i] == 53)
            {
                localU[i] = 0;
                localV[i] = 0;
            }
        }*/

        localMatrix.convectiveMembers(localMatrix3, localU, localV, a, b);
        //localMatrix.supg.convectiveMembers(a, b, localU, localV, hElem, mesh.square[t], mu, localMatrix3);

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += localVector0[i] - (9.81 * mesh.square[t]) / 3.0;
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
        methods.multMC(localMatrix0, 2.0, 3);

        //localMatrix.derivative(localMatrix2, b);
        //methods.multMV(localMatrix2, localU, localVector0, 3);
        //methods.actionsVC(localVector0, 1.0 / del_t, 3, '*');

        localMatrix.derivative(localMatrix2, a);
        methods.multMV(localMatrix2, localV, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

        for (unsigned int i = 0; i < 3; i++)
        {
            B[k[i]] += 0;//localVector1[i]; //
            for (unsigned int j = 0; j < 3; j++)
            {
                A[k[i]][k[j]] += -localMatrix0[i][j];// - localMatrix1[i][j];// + localMatrix1[i][j];
            }
        }
    }

    /*for (int i = 0; i < mesh.n; i++)
    {
        double yy = mesh.nodes[i][0];

        if (yy == 0)
        {
            map<unsigned int, double>::iterator iterMap;
            for (iterMap = A[i].begin(); iterMap != A[i].end(); iterMap++)
            {
                //Tn[i] += (*iterMap).second * Sn[(*iterMap).first];
                cout << i+1 << ":" << (*iterMap).second << " ";
            }
            cout << endl;
        }

    }*/


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

void NavierStokes::fillSLAE_c()
{
    double hElem;
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

        localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], d, localMatrix0);
        methods.multMC(localMatrix0, 1.0 / del_t, 3);

        localMatrix.supg.convectiveMembers(a, b, localU, localV, hElem, mesh.square[t], d, localMatrix1);

        localMatrix.supg.mass(a, b, localU, localV, hElem, mesh.square[t], d, localMatrix2);
        methods.multMV(localMatrix2, localVector0, localVector1, 3);
        methods.actionsVC(localVector1, 1.0 / del_t, 3, '*');

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

double NavierStokes::calcH(double *u, double *v, double square)
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

void NavierStokes::bc45()
{
    int numberBorder = 3;
    double f;  // среднее значение функции на отрезке
    double X1, Y1, X2, Y2;   //координаты отрезка
    double l; //длинна отрезка
    unsigned int countBorderElements; //количество отрезков на данной границе
    unsigned int k1, k2;   //номера точек отрезка

    for (unsigned int k = 0; k < mesh.border[numberBorder].m; k++) //цикл по отрезку// {
    {
        k1 = mesh.border[numberBorder].line[k][0];     //
        k2 = mesh.border[numberBorder].line[k][1]; //  находим номера точек первого отрезка

        X1 = mesh.nodes[k1][0];  //
        Y1 = mesh.nodes[k1][1];  // берём координаты первого отрезка
        X2 = mesh.nodes[k2][0]; //
        Y2 = mesh.nodes[k2][1]; //

        double c1 = rho / (del_t * sqrt(2)) * (uStar[k1] - vStar[k1]);
        double c2 = rho / (del_t * sqrt(2)) * (uStar[k2] - vStar[k2]);

        l = methods.findLength(X1, X2, Y1, Y2); //нахождение длинны отрезка
        f = (c1 + c2) / 2; //нахождение среднего значения функции на отрезке
        B[k1] -= l / 2 * f;
        B[k2] -= l / 2 * f;
    }

}




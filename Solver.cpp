//
// Created by konstantin on 27.02.19.
//

#include "Solver.h"

Solver::Solver(Mesh mesh, Arguments arguments, string pathToResult)
{
    this->mesh = mesh;
    this->arguments = arguments;
    this->pathToResult = pathToResult;

    reservMemory();
    init();
    setSquare();
}

void Solver::reservMemory()
{
    k = new int[3];
    x = new double[3];
    y = new double[3];
    a = new double[3];
    b = new double[3];

    A.resize(mesh.n);
    B = new double[mesh.n];

    mesh.square = new double[mesh.m];
}

void Solver::init()
{
    methods.null(B, mesh.n);

    h = mesh.h;

    startTime = arguments.time.start;
    endTime = arguments.time.end;

    if (arguments.time.del_t == 0)
        del_t = h - h / 5;
    else
        del_t = arguments.time.del_t;

    intervalRecord = arguments.write.interval;

    runTime = 0;
    tempTime = 1;
    countIterations = 1;
}

void Solver::setSquare()
{
    for (unsigned int t = 0; t < mesh.m; t++)
    {
        for (unsigned int i = 0; i < 3; i++)
        {
            k[i] = mesh.elements[t][i];
        }

        for (unsigned int i = 0; i<3; i++)
        {
            x[i] = mesh.nodes[k[i]][0];
            y[i] = mesh.nodes[k[i]][1];
        }

        a[0] = x[2] - x[1];   a[1] = x[0] - x[2];   a[2] = x[1] - x[0];
        b[0] = y[1] - y[2];   b[1] = y[2] - y[0];   b[2] = y[0] - y[1];

        mesh.square[t] = (b[0] * a[1] - a[0] * b[1]) / 2.0;
    }
}

void Solver::solveSLAE(double *vectorResult)
{
    methods.BISC_Mod(A, B, vectorResult, mesh.n);
    for (unsigned int i = 0; i < mesh.n; i++)
    {
        A[i].clear();
    }
    methods.null(B, mesh.n);
}

void Solver::autoRecordData(vector<double *> column, vector<string> name)
{

    if (arguments.write.record)
    {
        if (!arguments.write.onlyLastFrame)
        {
            if (runTime - intervalRecord * tempTime >= 0)
            {
                methods.recordDataTecplot(
                        mesh.nodes,
                        mesh.n,
                        mesh.elements,
                        mesh.m,
                        column,
                        name,
                        to_string(runTime),
                        pathToResult,
                        runTime
                );

                tempTime++;
            }


        } else
        {
            if (runTime >= endTime)
            {
                methods.recordDataTecplot(
                        mesh.nodes,
                        mesh.n,
                        mesh.elements,
                        mesh.m,
                        column,
                        name,
                        to_string(runTime) + "_last",
                        pathToResult,
                        runTime
                        );
            }
        }
    }
}

void Solver::recordData(vector<double *> column, vector<string> name)
{
    if (arguments.write.record && !arguments.write.onlyLastFrame)
    {
        methods.recordDataTecplot(
                mesh.nodes,
                mesh.n,
                mesh.elements,
                mesh.m,
                column,
                name,
                to_string(runTime),
                pathToResult,
                runTime
        );
    }
}

void Solver::conditionNode_1(double c, unsigned int dot)
{
    for (unsigned int i = 0; i < mesh.n; i++)   //цикл по вектору правой часки
    {
        B[i] -= c * A[i][dot];
    }

    B[dot] = c;  // подстановка граничных условий

    A[dot].clear();
    for (unsigned int i = 0; i < mesh.n; i++)
    {
        A[i].erase(dot);
    }
    A[dot][dot] = 1;
}

void Solver::conditionBorder_1(double c, unsigned int numberBorder)
{
    unsigned int t; //точка на границе
    unsigned int countDotBorder; //количество точек на границе
    for (unsigned int i = 0; i < mesh.border[numberBorder].n; i++)  //цикл по границе
    {
        t = mesh.border[numberBorder].numberNodes[i]; //выбираем точку на границе
        conditionNode_1(c, t);
    }
}

void Solver::conditionBorder_2(double c, unsigned int numberBorder)
{
    double f;  // среднее значение функции на отрезке
    double X1, Y1, X2, Y2;   //координаты отрезка
    double l; //длинна отрезка
    unsigned int countBorderElements; //количество отрезков на данной границе
    unsigned int k1, k2;   //номера точек отрезка

    for (unsigned int k = 0; k < mesh.border[numberBorder].m; k++) //цикл по отрезку
    {
        k1 = mesh.border[numberBorder].line[k][0];     //
        k2 = mesh.border[numberBorder].line[k][1]; //  находим номера точек первого отрезка

        X1 = mesh.nodes[k1][0];  //
        Y1 = mesh.nodes[k1][1];  // берём координаты первого отрезка
        X2 = mesh.nodes[k2][0]; //
        Y2 = mesh.nodes[k2][1]; //
        l = methods.findLength(X1, X2, Y1, Y2); //нахождение длинны отрезка
        f = (c + c) / 2; //нахождение среднего значения функции на отрезке
        B[k1] -= l / 2 * f;
        B[k2] -= l / 2 * f;
    }
}

void Solver::timeCutdown()
{
    if (runTime == 0)
    {
        cout << "0.0 %" << endl;
    } else
    {
        cout << (endTime / runTime) * 100 << " %" << endl;
    }

}

void Solver::calc_a_b()
{
    a[0] = x[2] - x[1];   a[1] = x[0] - x[2];   a[2] = x[1] - x[0];
    b[0] = y[1] - y[2];   b[1] = y[2] - y[0];   b[2] = y[0] - y[1];
}

void Solver::dc_dx_zero(double **localcMatrix, int numberBorder)
{
    /*int t = mesh.border[numberBorder].numberNodes[0]; //select nodes on border
    for (int i = 0; i < 3; i++)
    {
        if (x[i] == mesh.nodes[t][0])
        {
            localcMatrix[i][0] = 0;
            localcMatrix[i][1] = 0;
            localcMatrix[i][2] = 0;
        }
    }*/
}

void Solver::dc_dy_zero(double **localcMatrix, int numberBorder)
{
    /*int t = mesh.border[numberBorder].numberNodes[0]; //select nodes on border
    for (int i = 0; i < 3; i++)
    {
        if (y[i] == mesh.nodes[t][1])
        {
            localcMatrix[i][0] = 0;
            localcMatrix[i][1] = 0;
            localcMatrix[i][2] = 0;
        }
    }*/
}



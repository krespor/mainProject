//
// Created by konstantin on 27.02.19.
//

#ifndef NEWFEM_SOLVER_H
#define NEWFEM_SOLVER_H

#include "./structurs/mesh.h"
#include "./structurs/arguments.h"
#include "./structurs/localMatrix.h"
#include "./structurs/methods.h"

#include <iostream>
#include "math.h"

using namespace std;



class Solver {

public:
    Solver(Mesh mesh, Arguments arguments, string pathToResult);

protected:
    Methods methods;

    Mesh mesh;
    Arguments arguments;
    LocalMatrix localMatrix;

    string pathToResult;

    double *x, *y, *a, *b;
    int *k;

    int countIterations;

    float intervalRecord;
    float tempTime;
    float startTime, endTime, del_t, runTime;
    float h;

    vector<map<unsigned int, double>> A;
    double *B;

    void solveSLAE(double *vectorResult);
    void autoRecordData(vector<double*> column, vector<string> name);
    void recordData(vector<double*> column, vector<string> name);

    void calc_a_b();

    void conditionNode_1(double c, unsigned int dot);
    void conditionBorder_1(double c, unsigned int numberBorder);
    void conditionBorder_2(double c, unsigned int numberBorder);

    void dc_dx_zero(double **localcMatrix, int numberBorder);
    void dc_dy_zero(double **localcMatrix, int numberBorder);

private:
    void init();
    void reservMemory();
    void setSquare();
};

#endif //NEWFEM_SOLVER_H

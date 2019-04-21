//
// Created by konstantin on 28.02.19.
//

#ifndef NEWFEM_SUPG_H
#define NEWFEM_SUPG_H

#include "./Solver.h"

class SUPG : Solver{

public:
    SUPG(Mesh mesh, Arguments arguments, string pathToResult);

private:

    double *c, *cn, *u, *v;
    double **localMatrix0, **localMatrix1, **localMatrix2, **localMatrix3;
    double *localVector0, *localVector1, *localVector2, *localVector3;
    double *localU, *localV;

    double hElem, d;

    void reservMemory();
    void init();

    void fillSLAE();

    void calc();
    double calcH(double *u, double *v, double square);
    double calcH_old(double *u, double *v, double square);

    void supgMatrixMass(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix);
    void supgMatrixLaplas(double *u, double *v, double h, double square, double k, double **matrix);
    void supgFull(double *a, double *b, double *u, double *v, double h, double square, double k, double **matrix);


};


#endif //NEWFEM_SUPG_H

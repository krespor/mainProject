//
// Created by konstantin on 06.04.19.
//

#ifndef NEWFEM_TWOPHASEFLOW_H
#define NEWFEM_TWOPHASEFLOW_H

#include "./Solver.h"

class TwoPhaseFlow : Solver{

public:
    TwoPhaseFlow(Mesh mesh, Arguments arguments, string pathToResult);

private:
    double mu0, mu1, rho0, rho1, g;
    double *p, *u, *v, *un, *vn, *uStar, *vStar, *alpha, *alpha_n, *mu, *rho;

    double **localMatrix0, **localMatrix1, **localMatrix2, **localMatrix3;
    double *localVector0, *localVector1, *localVector2, *localVector3;
    double *localU, *localV, *localP, *localRho, *localMu, *localAlpha_n;

    void fillSLAE_uStar();
    void fillSLAE_testVstar();
    void fillSLAE_vStar();
    void fillSLAE_p();
    void fillSLAE_testP();
    void fillSLAE_testPavg();
    void fillSLAE_u();
    void fillSLAE_v();
    void fillSLAE_testV();
    void fillSLAE_alpha();

    void bcPressure();

    void currectUV();
    double calcH(double *u, double *v, double square);
    void setPhase();

    void reservMemory();
    void init();
    void calc();

};


#endif //NEWFEM_TWOPHASEFLOW_H

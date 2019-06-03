//
// Created by konstantin on 06.03.19.
//

#ifndef NEWFEM_NAVIERSTOKES_H
#define NEWFEM_NAVIERSTOKES_H

#include "./Solver.h"

class NavierStokes : Solver{
public:
    NavierStokes(Mesh mesh, Arguments arguments, string pathToResult);

private:

    double rho, mu, d;
    double *p, *u, *v, *un, *vn, *uStar, *vStar, *c, *cn;

    double **localMatrix0, **localMatrix1, **localMatrix2, **localMatrix3;
    double *localVector0, *localVector1, *localVector2, *localVector3;
    double *localU, *localV;

    void fillSLAE_uStar();
    void fillSLAE_vStar();
    void fillSLAE_p();
    void fillSLAE_u();
    void fillSLAE_v();
    void fillSLAE_c();

    void currectUV();
    void privateBorderCondition();
    double calcH(double *u, double *v, double square);

    void reservMemory();
    void init();
    void calc();

    void calcBackwardStep();
    void calcCaverna();
    void calcTwoPhaseFlow();

    void bc45();

};


#endif //NEWFEM_NAVIERSTOKES_H

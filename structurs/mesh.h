//
// Created by konstantin on 26.02.19.
//

#ifndef NEWFEM_MESH_H
#define NEWFEM_MESH_H

#include <iostream>

struct Border
{
    unsigned int n, m;
    unsigned int *numberNodes;
    unsigned int **line;
};

struct Mesh
{
    std::string name;
    float h;
    unsigned int n, m;
    double **nodes;
    unsigned  int **elements;
    double *square;

    int countBorder;
    Border *border;
};

#endif //NEWFEM_MESH_H

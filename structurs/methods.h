//
// Created by konstantin on 27.02.19.
//

#ifndef NEWFEM_METHODS_H
#define NEWFEM_METHODS_H

#include "math.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>

using namespace std;

struct Methods
{
    void multMV(double **a, double *b, double *c, int n)
    {
        for (int i = 0; i<n; i++) {
            c[i] = 0;
            for (int j = 0; j<n; j++) {
                c[i] += (a[i][j] * b[j]);
            }
        }
    }

    void multMC(double ** a, double c, int n)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a[i][j] *= c;
            }
        }
    }

    void null(double *a, int n)
    {
        for (int i = 0; i<n; i++)
            a[i] = 0;
    }

    double maxElement(double * a, int n)
    {
        double max = a[0];
        for (int i = 1; i < n; i++)
        {
            if (a[i] > max)
            {
                max = a[i];
            }
        }
        return max;
    }

    double minElement(double * a, int n)
    {
        double min = a[0];
        for (int i = 1; i < n; i++)
        {
            if (a[i] < min)
            {
                min = a[i];
            }
        }
        return min;
    }

    double calcMistake(double * a, double * b, int n)
    {
        double *error = new double[n];
        double tmp;

        for (int i = 0; i < n; i++)
        {
            error[i] = fabs(a[i] - b[i]);
        }

        tmp = maxElement(error, n);

        delete[]error;

        return tmp;
    }

    double scalar(double * a, double * b, int n)
    {
        double s = 0;
        for (int i = 0; i < n; i++)
        {
            s += a[i] * b[i];
        }
        return s;
    }

    void equateV(double * a, double * b, int n)
    {
        for (int i = 0; i < n; i++)
        {
            a[i] = b[i];
        }
    }

    void null(double **a, int n)
    {
        for (int i = 0; i<n; i++)
            for (int j = 0; j<n; j++)
                a[i][j] = 0;
    }

    double findLength(double x1, double x2, double y1, double y2)
    {
        return sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
    }

    void actionsVC(double * a, double c, int n, char act)
    {
        switch (act)
        {
            case '+':
                for (int i = 0; i < n; i++)
                {
                    a[i] += c;
                }
                break;
            case '-':
                for (int i = 0; i < n; i++)
                {
                    a[i] -= c;
                }
                break;
            case '*':
                for (int i = 0; i < n; i++)
                {
                    a[i] *= c;
                }
                break;
            case '/':
                for (int i = 0; i < n; i++)
                {
                    a[i] /= c;
                }
                break;
            case '=':
                for (int i = 0; i < n; i++)
                {
                    a[i] = c;
                }
                break;
            default:
                break;
        }
    }

    int BISC_Mod(vector<map<unsigned int, double>> AA, double * B, double * X, int n)
    {
        map<unsigned int, double>::iterator iterMap;

        double *Rn, *_Rn, *Pn, *Pn1, *Vn, *A, eps = 10e-10, *Sn, *Tn;
        int iter = 0;
        Rn = new double[n];
        Pn = new double[n];
        Pn1 = new double[n];
        _Rn = new double[n];
        Vn = new double[n];
        A = new double[n];
        Sn = new double[n];
        Tn = new double[n];

        null(A, n);

        for (int i = 0; i < n; i++)
        {
            X[i] = B[i] + 0.000001;
        }


        for (int i = 0; i < n; i++)
        {
            for (iterMap = AA[i].begin(); iterMap != AA[i].end(); iterMap++)
            {
                A[i] += (*iterMap).second * X[(*iterMap).first];  //������������ � �� 1�� ����������� ������� �=1
            }
            _Rn[i] = Rn[i] = B[i] - A[i];
        }

        null(Pn, n);
        null(Vn, n);

        double pn = 1.;
        double an = 1.;
        double bn = 1.;
        double wn = 1.;
        double pn1;
        double norm = 0;

        int p = 1;
        while (true)
        {
            pn1 = scalar(_Rn, Rn, n);
            //
            if (pn1 == 0)
            {
                //cout << "ERROR!" << endl;
            }
            if (p == 1)
            {
                for (int i = 0; i < n; i++)
                {
                    Pn1[i] = Rn[i];
                }
            }
            else
            {
                bn = (pn1 / pn)*(an / wn);
                for (int i = 0; i < n; i++)
                {
                    Pn1[i] = Rn[i] + bn * (Pn[i] - wn * Vn[i]);
                }
            }
            //
            null(Vn, n);
            for (int i = 0; i < n; i++)
            {
                for (iterMap = AA[i].begin(); iterMap != AA[i].end(); iterMap++)
                {
                    Vn[i] += (*iterMap).second * Pn1[(*iterMap).first];
                }
            }

            an = pn1 / scalar(_Rn, Vn, n);
            for (int i = 0; i < n; i++)
            {
                Sn[i] = Rn[i] - an * Vn[i];
            }
            //cout << "n == " << sqrt(scalar(Sn, Sn, n)) << endl;
            if ((sqrt(scalar(Sn, Sn, n))) <= eps)
            {
                for (int i = 0; i < n; i++)
                {
                    X[i] += an * Pn1[i];
                }

                delete[]Rn;
                delete[]Pn;
                delete[]Pn1;
                delete[]_Rn;
                delete[]Vn;
                delete[]A;
                delete[]Sn;
                delete[]Tn;

                return 1;
            }
            else
            {
                null(Tn, n);
                for (int i = 0; i < n; i++)
                {
                    for (iterMap = AA[i].begin(); iterMap != AA[i].end(); iterMap++)
                    {
                        Tn[i] += (*iterMap).second * Sn[(*iterMap).first];
                    }
                }

                wn = scalar(Tn, Sn, n) / scalar(Tn, Tn, n);
                for (int i = 0; i < n; i++)
                {
                    X[i] += an * Pn1[i] + wn * Sn[i];
                    Rn[i] = Sn[i] - wn * Tn[i];
                }

                pn = pn1;
                for (int i = 0; i < n; i++)
                {
                    Pn[i] = Pn1[i];
                }
                //cout << "n = " << sqrt(scalar(Rn, Rn, n)) << endl;
                if ((sqrt(scalar(Rn, Rn, n))) <= eps)
                {

                    delete[]Rn;
                    delete[]Pn;
                    delete[]Pn1;
                    delete[]_Rn;
                    delete[]Vn;
                    delete[]A;
                    delete[]Sn;
                    delete[]Tn;

                    return 1;
                }
            }
            p++;
            iter++;
        }
    }
};

#endif //NEWFEM_METHODS_H

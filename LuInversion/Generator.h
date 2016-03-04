#pragma once

#include <math.h>
#include <iostream>
#include "Matrix.h"

#define SIGN_LAW 1
#define LAMBDA_LAW 2
#define VARIANT 1
#define SCHEMA 1

using namespace std;

class Generator {
private:
    void Q_matrix(double **Q, int n, int schema);

    void matr_mul(double **a, double **b, double **c, int n);

    double matr_inf_norm(double **a, int n);

    void mygen(double **a, double **a_inv, int n, double alpha, double beta, int sign_law, int lambda_law, int variant,
               int schema);

    Matrix a;
    Matrix aInv;
    double norm_inv;
    double norm;
    double obusl;
    double Rnorm;

public:
    Generator(int n, double alpha, double beta);

    Matrix getMatrix();

    Matrix getInvertedMatrix();

    double getNorm_inv() const;

    double getNorm() const;

    double getObusl() const;

    double getRnorm() const;
};
#pragma once

#include <iostream>
#include <string.h>
#include "Matrix.h"
#include <math.h>

//using namespace std;

class Vector {
private:
    double *vector;
    size_t size;
public:
    Vector(size_t N);

    Vector(const Vector &);

    Vector &operator=(const Vector &);

    ~Vector();

    size_t getSize() const;

    Vector operator-(const Vector &) const;

    double Norma() const;

    double &operator()(size_t i);

    double operator()(size_t i) const;

    friend std::ostream &operator<<(std::ostream &os, const Vector &);

    friend Vector operator*(double, const Vector &);
};

Vector operator*(const Matrix &, const Vector &);
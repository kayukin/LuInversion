#pragma once

#include <stddef.h>
#include <iosfwd>
#include <string.h>
#include <ostream>
#include <cmath>
#include <istream>
using namespace std;

class Matrix {
protected:
    double **matr;
    size_t size;

    void Destroy();

    void Copy(const Matrix &);

public:
    Matrix(size_t N);

    Matrix(const Matrix &);

    Matrix(istream& ifs);

    ~Matrix();

    Matrix &operator=(const Matrix &);

    Matrix operator*(const Matrix &B) const;

	Matrix operator-(const Matrix &B)const;

    double Norm() const;

    double Determinant() const;

    friend ostream &operator<<(ostream &os, const Matrix &matrix);

    double &operator()(size_t i, size_t j);

    size_t getSize() const;

    double operator()(size_t i, size_t j) const;

    friend class LuDecomposition;
	friend class Generator;
};
#ifndef MINERRMETHOD_MATRIX_H
#define MINERRMETHOD_MATRIX_H

#include <stddef.h>
#include <iosfwd>
#include <string.h>
#include <ostream>
#include <cmath>

class Matrix {
protected:
    double **matr;
    size_t size;

    void Destroy();

    void Copy(const Matrix &);

public:
    Matrix(size_t N);

    Matrix(const Matrix &);

    ~Matrix();

    Matrix &operator=(const Matrix &);

    Matrix operator*(const Matrix &B) const;

    double Norm() const;

    double Determinant() const;

    friend std::ostream &operator<<(std::ostream &os, const Matrix &matrix);

    double &operator()(size_t i, size_t j);

    size_t getSize() const;

    double operator()(size_t i, size_t j) const;

    friend class LuDecomposition;
};

#endif //MINERRMETHOD_MATRIX_H

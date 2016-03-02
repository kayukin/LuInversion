//
// Created by Konstantin on 29.02.2016.
//

#include <iostream>
#include "LuDecomposition.h"

void PII_LUDecomposition(double **, int);

Matrix LuDecomposition::Inverse(const Matrix &matrix) {
    size_t size = matrix.getSize();
    Matrix invertedMatrix(size);
    Matrix LU(size);
    Matrix LUPII(matrix);
    PII_LUDecomposition(LUPII.matr, size);

    double **luPointer = LU.matr;
    double **invertedMatrixPointer = invertedMatrix.matr;
    double **matrixPointer = matrix.matr;


    //выбор максимального по модулю
    //использовать одну матрицу

    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            double sum = 0;
            if (i <= j) {
                for (size_t k = 0; k < i; k++)
                    sum += luPointer[i][k] * luPointer[k][j];
                luPointer[i][j] = matrixPointer[i][j] - sum;
            }
            else {
                for (size_t k = 0; k < j; k++)
                    sum += luPointer[i][k] * luPointer[k][j];
                if (luPointer[j][j] == 0)
                    return Matrix(0);
                luPointer[i][j] = (matrixPointer[i][j] - sum) / luPointer[j][j];
            }
        }
    }

    std::cout << LU << std::endl;
    std::cout << LUPII << std::endl;

    for (int i = size - 1; i >= 0; i--) {
        for (int j = size - 1; j >= 0; j--) {
            double sum = 0;
            if (i == j) {
                for (int p = j + 1; p < size; p++)
                    sum += luPointer[j][p] * invertedMatrixPointer[p][j];
                invertedMatrixPointer[j][j] = (1 - sum) / luPointer[j][j];
            }
            else if (i < j) {
                for (size_t p = i + 1; p < size; p++)
                    sum += luPointer[i][p] * invertedMatrixPointer[p][j];
                invertedMatrixPointer[i][j] = -sum / luPointer[i][i];
            }
            else {
                for (int p = j + 1; p < size; p++)
                    sum += invertedMatrixPointer[i][p] * luPointer[p][j];
                invertedMatrixPointer[i][j] = -sum;
            }

        }
    }

    return invertedMatrix;
}

void PII_LUDecomposition(double **m, int n) {
/*
Factors "m" matrix into A=LU where L is lower triangular and U is upper
triangular. The matrix is overwritten by LU with the diagonal elements
of L (which are unity) not stored. This must be a square n x n matrix.
*/

    int i;
    int *ri = new int[n];

    // Initialize the pointer vector.

    for (i = 0; i < n; i++)
        ri[i] = i;

    // LU factorization.

    for (int p = 1; p <= n - 1; p++) {
        // Find pivot element.
        for (i = p + 1; i <= n; i++) {
            if (fabs(m[ri[i - 1]][p - 1]) > fabs(m[ri[p - 1]][p - 1])) {
                // Switch the index for the p-1 pivot row if necessary.
                int t = ri[p - 1];
                ri[p - 1] = ri[i - 1];
                ri[i - 1] = t;
            }
        }

        if (m[ri[p - 1]][p - 1] == 0) {
            // The matrix is singular.
            return;
        }

        // Multiply the diagonal elements.
        // Form multiplier.

        for (i = p + 1; i <= n; i++) {
            m[ri[i - 1]][p - 1] /= m[ri[p - 1]][p - 1];

            // Eliminate [p-1].

            for (int j = p + 1; j <= n; j++)
                m[ri[i - 1]][j - 1] -= m[ri[i - 1]][p - 1] * m[ri[p - 1]][j - 1];

        }

    }
    delete[] ri;
}
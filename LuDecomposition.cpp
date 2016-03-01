//
// Created by Konstantin on 29.02.2016.
//

#include <iostream>
#include "LuDecomposition.h"

Matrix LuDecomposition::Inverse(const Matrix &matrix) {
    size_t size = matrix.getSize();
    Matrix invertedMatrix(size);
    Matrix LU(size);

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
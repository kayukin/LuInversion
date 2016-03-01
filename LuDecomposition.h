//
// Created by Konstantin on 29.02.2016.
//

#ifndef LU_INVERSION_LUDECOMPOSITION_H
#define LU_INVERSION_LUDECOMPOSITION_H


#include "Matrix.h"

class LuDecomposition {

public:
    static Matrix Inverse(const Matrix &matrix);
};


#endif //LU_INVERSION_LUDECOMPOSITION_H

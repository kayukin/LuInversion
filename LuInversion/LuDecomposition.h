#pragma once

#include "Matrix.h"

class LuDecomposition {

public:
    static Matrix Inverse(const Matrix &matrix);
private:
	static size_t findMaxColumnIndex(double* row, size_t size);
	static void swapColumns(double** matrix, size_t size, size_t first, size_t second, size_t startRow);
};
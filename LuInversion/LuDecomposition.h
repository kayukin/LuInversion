#pragma once
#include "Matrix.h"

class LuDecomposition
{
public:
	static Matrix Inverse(Matrix& matrix);
private:
	static void luDecompose(Matrix& matrix);
};

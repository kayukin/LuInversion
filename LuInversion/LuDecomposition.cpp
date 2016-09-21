#include "LuDecomposition.h"


struct DecompositionExeption
{
};

Matrix LuDecomposition::Inverse(Matrix& matrix)
{
	size_t size = matrix.getSize();
	Matrix invertedMatrix(size);

	luDecompose(matrix);

	//LU = matrix;
	//luDecomposeSingleMatrix(luPointer, size);


	for (int i = size - 1; i >= 0; i--)
	{
		for (int j = size - 1; j >= 0; j--)
		{
			double sum = 0;
			if (i == j)
			{
				for (int p = j + 1; p < size; p++)
					sum += matrix(j, p) * invertedMatrix(p, j);
				invertedMatrix(j, j) = (1 - sum) / matrix(j, j);
			}
			else if (i < j)
			{
				for (size_t p = i + 1; p < size; p++)
					sum += matrix(i, p) * invertedMatrix(p, j);
				invertedMatrix(i, j) = -sum / matrix(i, i);
			}
			else
			{
				for (int p = j + 1; p < size; p++)
					sum += invertedMatrix(i, p) * matrix(p, j);
				invertedMatrix(i, j) = -sum;
			}
		}
	}

	return invertedMatrix;
}

void LuDecomposition::luDecompose(Matrix& matrix)
{
	for (size_t i = 0; i < matrix.getSize(); ++i)
	{
		for (size_t j = 0; j < matrix.getSize(); ++j)
		{
			double sum = 0;
			if (i <= j)
			{
				for (size_t k = 0; k < i; k++)
					sum += matrix(i, k) * matrix(k, j);
				matrix(i, j) = matrix(i, j) - sum;
			}
			else
			{
				for (size_t k = 0; k < j; k++)
					sum += matrix(i, k) * matrix(k, j);
				if (matrix(j, j) == 0)
					throw DecompositionExeption();
				matrix(i, j) = (matrix(i, j) - sum) / matrix(j, j);
			}
		}
	}
}

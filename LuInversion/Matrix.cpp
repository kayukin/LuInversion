#include "Matrix.h"

#include <fstream>

using namespace std;

Matrix::Matrix(size_t N)
{
	size = N;
	matr = new double *[size];
	for (size_t i = 0; i < size; i++)
	{
		matr[i] = new double[size];
		memset(matr[i], 0, sizeof(double) * size);
	}
}

Matrix::Matrix(const Matrix& matrix)
{
	Copy(matrix);
}

Matrix::~Matrix()
{
	Destroy();
}

Matrix& Matrix::operator=(const Matrix& matrix)
{
	if (this == &matrix)
	{
		return *this;
	}
	Destroy();
	Copy(matrix);
	return *this;
}


void Matrix::Destroy()
{
	for (size_t i = 0; i < size; i++)
		delete[] matr[i];
	delete[] matr;
}

void Matrix::Copy(const Matrix& matrix)
{
	size = matrix.size;
	matr = new double *[size];
	for (size_t i = 0; i < size; i++)
	{
		matr[i] = new double[size];
		for (size_t j = 0; j < size; j++)
		{
			matr[i][j] = matrix.matr[i][j];
		}
	}
}

Matrix Matrix::operator*(const Matrix& B) const
{
	Matrix result(size);
	for (size_t i = 0; i < size; i++)
		for (size_t j = 0; j < size; j++)
		{
			for (size_t k = 0; k < size; k++)
				result.matr[i][j] += this->matr[i][k] * B.matr[k][j];
		}
	return result;
}

Matrix Matrix::operator-(const Matrix& B) const
{
	Matrix res(*this);
	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			res.matr[i][j] -= B.matr[i][j];
		}
	}
	return res;
}

double Matrix::Norm() const
{
	double s, norm = 0.;
	for (size_t i = 0; i < size; i++)
	{
		s = 0.;
		for (size_t j = 0; j < size; j++)
			s += fabs(matr[i][j]);
		if (s > norm)
			norm = s;
	}
	return norm;
}

double Matrix::Determinant() const
{
	int l;
	double d;
	double sum11 = 1, sum12 = 0, sum21 = 1, sum22 = 0;
	for (int i = 0; i < size; i++)
	{
		sum11 = 1;
		l = (int) (2 * size - 1 - i);
		sum21 = 1;
		for (int j = 0; j < size; j++)
		{
			sum21 *= matr[j][l % size];
			l--;
			sum11 *= matr[j][(j + i) % (size)];
		}
		sum22 += sum21;
		sum12 += sum11;
	}
	d = sum12 - sum22;
	return d;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix)
{
	for (size_t i = 0; i < matrix.size; i++)
	{
		for (size_t j = 0; j < matrix.size; j++)
			os << matrix.matr[i][j] << '\t';
		os << std::endl;
	}
	return os;
}

size_t Matrix::getSize() const
{
	return size;
}

Matrix::Matrix(istream& ifs)
{
	ifs >> size;
	matr = new double *[size];
	for (size_t i = 0; i < size; i++)
	{
		matr[i] = new double[size];
		memset(matr[i], 0, sizeof(double) * size);
	}
	for (size_t i = 0; i < size; i++)
		for (size_t j = 0; j < size; j++)
		{
			ifs >> matr[i][j];
		}
}

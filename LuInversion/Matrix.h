#pragma once

#include <iosfwd>
#include <ostream>


class Matrix
{
protected:
	double** matr;
	size_t size;

	void Destroy();

	void Copy(const Matrix&);

public:
	explicit Matrix(size_t N);
	explicit Matrix(std::istream& ifs);

	Matrix(const Matrix& other);
	Matrix& operator=(const Matrix& other);

	~Matrix();


	Matrix operator*(const Matrix& B) const;

	Matrix operator-(const Matrix& B) const;

	double Norm() const;

	double Determinant() const;

	friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);

	inline double& operator()(size_t i, size_t j)
	{
		return matr[i][j];
	}

	inline double operator()(size_t i, size_t j) const
	{
		return matr[i][j];
	}

	size_t getSize() const;


	friend class LuDecomposition;
	friend class Generator;
};

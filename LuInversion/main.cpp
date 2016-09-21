#include <iostream>
#include <fstream>
#include "Generator.h"
#include "LuDecomposition.h"

#define N 3


using namespace std;

int main()
{
	double beta = 10;
	double alpha = 1;
	size_t n = N;

	ofstream ofs("output.csv");
	ofstream txt("output.txt");

	ofs << "alpha, beta, Norm from generator, Norm, obusl" << endl;
	for (beta = 10; beta < 1e4; beta *= 10)
	{
		cout << '.';
		Generator generator(n, alpha, beta);
		cout << '.';
		Matrix invertedMatrix = LuDecomposition::Inverse(generator.getMatrix());
		cout << '.';
		ofs << alpha << ',' << beta << ',' << generator.getNorm_inv()
			<< ',' << invertedMatrix.Norm() << ',' << generator.getObusl() << endl;
		txt << invertedMatrix << endl << endl;
		txt << generator.getInvertedMatrix() << endl << endl;
		txt << "-----------------------------------------------" << endl << endl;
	}
	beta = 10;
	for (alpha = 1e-1; alpha > 1e-4; alpha /= 10)
	{
		cout << '.';
		Generator generator(n, alpha, beta);
		cout << '.';
		Matrix invertedMatrix = LuDecomposition::Inverse(generator.getMatrix());
		cout << '.';
		ofs << alpha << ',' << beta << ',' << generator.getNorm_inv()
			<< ',' << invertedMatrix.Norm() << ',' << generator.getObusl() << endl;
	}
}

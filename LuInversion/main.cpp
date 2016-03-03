#include <iostream>
#include <fstream>
#include "Generator.h"
#include "LuDecomposition.h"

using namespace std;

int main() {
    double beta = 10;
    double alpha = 1;
    size_t n = 3;

    ofstream ofs("output.csv");
	ofstream txt("output.txt");

/*
	ofs << "Generated matrix" << endl;
	ofs << generator.getMatrix() << endl;
	ofs << "Inverted matrix" << endl;
	ofs << invertedMatrix << endl;
	ofs << "Generated inv matrix" << endl;
	ofs << generator.getInvertedMatrix() << endl;
*/
	ofs << "alpha, beta, Norm from generator, Norm, obusl" << endl;
	for (beta = 10; beta < 1e4; beta *= 10) {
		Generator generator(n, alpha, beta);
		Matrix invertedMatrix = LuDecomposition::Inverse(generator.getMatrix());
		ofs << alpha << ',' << beta << ',' << generator.getNorm_inv() 
			<< ',' << invertedMatrix.Norm() << ',' << generator.getObusl() << endl;
		txt << invertedMatrix << endl << endl;
		txt << generator.getInvertedMatrix() << endl << endl;
		txt << "-----------------------------------------------" << endl << endl;
	}
	beta = 10;
	for (alpha = 1e-1; alpha > 1e-4; alpha /= 10) {
		Generator generator(n, alpha, beta);
		Matrix invertedMatrix = LuDecomposition::Inverse(generator.getMatrix());
		ofs << alpha << ',' << beta << ',' << generator.getNorm_inv()
			<< ',' << invertedMatrix.Norm() << ',' << generator.getObusl() << endl;
	}
}
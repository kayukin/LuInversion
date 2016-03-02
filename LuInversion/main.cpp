#include <iostream>
#include <fstream>
#include "Generator.h"
#include "LuDecomposition.h"

using namespace std;

int main() {
    double beta = 10;
    double alpha = 1;
    size_t n = 3;

    Generator generator(n, alpha, beta);
    ofstream ofs("output.txt");
	ofs << "Generated matrix" << endl;
	ofs << generator.getMatrix << endl;
    Matrix invertedMatrix = LuDecomposition::Inverse(generator.getMatrix());
	ofs << "Inverted matrix" << endl;
	ofs << invertedMatrix << endl;
	ofs << "Generated inv matrix" << endl;
	ofs << generator.getInvertedMatrix << endl;
}
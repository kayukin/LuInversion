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
    ifstream input("input.txt");
    Matrix myMatrix(input);
    Matrix invertedMatrix = LuDecomposition::Inverse(generator.getMatrix());

    cout << invertedMatrix << endl;
    cout << generator.getInvertedMatrix() << endl;
}
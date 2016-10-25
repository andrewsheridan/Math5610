//Andrew Sheridan
//Math 5610 
//Written in C++
#include "../../HW4/BackSubstitution/BackSubstitution/Matrix.h"
#include <iostream>

int main(void) {
	const unsigned size = 4;

	//Problem 1
	double** matrix = CreateMatrix(size);
	double* vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	double** U = CopyMatrix(matrix, size);
	double** L = ScaledLUFactorization(U, vector, size);

	std::cout << "Problem 1: Scaled Partial Pivoting" << std::endl;
	std::cout << "Matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of Scaled LU Factorization " << std::endl;
	PrintMatrix(U, size);
	PrintMatrix(L, size);

	int input;
	std::cin >> input;

	return 0;
}

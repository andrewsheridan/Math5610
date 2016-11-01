#include "Matrix.h";
#include "Vector.h"
#include <iostream>

int main() 
{
	const int size = 10;
	double** matrix = CreateSymmetricMatrix(size);

	if (IsMatrixSymmetric(matrix, size)) {
		std::cout << "Symmetric! 1" << std::endl;
	}

	matrix = CreateIdentityMatrix(size);
	if (IsMatrixSymmetric(matrix, size)) {
		std::cout << "Symmetric! 2" << std::endl;
	}
	int input; 
	std::cin >> input;
	return 0;
}
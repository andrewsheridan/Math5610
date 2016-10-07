//Andrew Sheridan
//Math 5610 
//Written in C++

//Main.cpp
#include<iostream>

#include "Matrix.h"


int main(void) {
	const unsigned size = 3;

	//Problem 1
	double** matrix = CreateUpperTriangularMatrix(size);
	double* vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	double** matrixCopy = CopyMatrix(matrix, size);
	double* result = BackSubstitution(matrixCopy, vector, size);

	std::cout << "Problem 1" << std::endl;
	std::cout << "Upper triangular matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of back substition " << std::endl;
	PrintVector(result, size);


	//Problem 2
	matrix = CreateLowerTriangularMatrix(size);
	vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	matrixCopy = CopyMatrix(matrix, size);
	result = ForwardSubstitution(matrixCopy, vector, size);

	std::cout << "Problem 2" << std::endl;
	std::cout << "Lower triangular matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of forward substition " << std::endl;
	PrintVector(result, size);

	//Problem 5
	
	matrix = CreateMatrix(size);
	vector = new double[size];
	vector[0] = 1;
	vector[1] = 2;
	vector[2] = 3;


	double* fivePartOne = VectorMatrixMultiply(matrix, vector, size);
	std::cout << "Problem 5" << std::endl;
	std::cout << "5.1" << std::endl;
	PrintMatrix(matrix, size);
	PrintVector(vector, size);
	PrintVector(fivePartOne, size);


	double* onesVector = CreateOnesVector(size);
	double* fivePartThree = VectorMatrixMultiply(matrix, onesVector, size);

	std::cout << "Results from 5.3" << std::endl;
	PrintVector(onesVector, size);
	PrintVector(fivePartThree, size);

	//Stop console from closing
	int input;
	std::cin >> input;

	return 0;
}
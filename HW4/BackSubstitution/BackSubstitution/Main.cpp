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
	double* resultVector = BackSubstitution(matrixCopy, vector, size);

	std::cout << "Problem 1: Back Substitution" << std::endl;
	std::cout << "Upper triangular matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of back substition " << std::endl;
	PrintVector(resultVector, size);


	//Problem 2
	matrix = CreateLowerTriangularMatrix(size);
	vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	matrixCopy = CopyMatrix(matrix, size);
	resultVector = ForwardSubstitution(matrixCopy, vector, size);

	std::cout << "Problem 2: Forward Substitution" << std::endl;
	std::cout << "Lower triangular matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of forward substition " << std::endl;
	PrintVector(resultVector, size);

	//Problem 3
	matrix = CreateDiagonallyDominantMatrix(size);
	vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	matrixCopy = CopyMatrix(matrix, size);
	resultVector = CopyVector(vector, size);
	double** resultMatrix = GaussianElimination(matrixCopy, resultVector, size);

	std::cout << "Problem 3: Gaussian Elimination" << std::endl;
	std::cout << "Test Matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of gaussian elimination " << std::endl;
	PrintAugmentedMatrix(resultMatrix, resultVector, size);

	//Problem 4
	matrix = CreateDiagonallyDominantMatrix(size);
	vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	matrixCopy = CopyMatrix(matrix, size);
	resultVector = CopyVector(vector, size);
	resultMatrix = GaussianElimination(matrixCopy, resultVector, size);

	double* finalVector = BackSubstitution(resultMatrix, resultVector, size);

	std::cout << "Problem 4: Gaussian Elimination w/ Back Substitution" << std::endl;
	std::cout << "Test Matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of gaussian elimination and back substitution" << std::endl;
	PrintAugmentedMatrix(resultMatrix, finalVector, size);

	//Problem 5
	matrix = CreateDiagonallyDominantMatrix(size);
	vector = new double[size];
	vector[0] = 1;
	vector[1] = 2;
	vector[2] = 3;


	double* fivePartOne = VectorMatrixMultiply(matrix, vector, size);
	std::cout << "Problem 5: Testing" << std::endl;
	std::cout << "--------------5.1--------------" << std::endl;
	PrintMatrix(matrix, size);
	PrintVector(vector, size);
	PrintVector(fivePartOne, size);


	double* onesVector = CreateOnesVector(size);
	double* fivePartThree = VectorMatrixMultiply(matrix, onesVector, size);

	std::cout << "--------------5.3--------------" << std::endl;
	PrintVector(onesVector, size);
	PrintVector(fivePartThree, size);

	std::cout << "--------------5.4--------------" << std::endl;	
	for (int i = 10; i *= 2; i <= 160) {
		matrix = CreateDiagonallyDominantMatrix(i);
		vector = CreateOnesVector(i);
		vector = VectorMatrixMultiply(matrix, vector, i);
		matrixCopy = CopyMatrix(matrix, i);
		resultMatrix = GaussianElimination(matrixCopy, vector, i);
		resultVector = BackSubstitution(resultMatrix, vector, i);

		/*std::cout << "Problem 5: Testing" << std::endl;
		std::cout << "Start matrix" << std::endl;
		PrintMatrix(matrix, i);
		std::cout << "Test vector " << std::endl;
		PrintVector(vector, i);*/
		std::cout << "Result of GE and Back Substitution " << std::endl;
		PrintVector(resultVector, i);
	}

	//Problem 4
	matrix = CreateDiagonallyDominantMatrix(size);
	double** matrixTwo = CreateIdentityMatrix(size);

	vector = CreateOnesVector(size);
	vector = VectorMatrixMultiply(matrix, vector, size);
	matrixCopy = CopyMatrix(matrix, size);
	resultVector = CopyVector(vector, size);
	resultMatrix = GaussianElimination(matrixCopy, resultVector, size);

	double* finalVector = BackSubstitution(resultMatrix, resultVector, size);

	std::cout << "Problem 4: Gaussian Elimination w/ Back Substitution" << std::endl;
	std::cout << "Test Matrix" << std::endl;
	PrintMatrix(matrix, size);
	std::cout << "Test vector " << std::endl;
	PrintVector(vector, size);
	std::cout << "Result of gaussian elimination and back substitution" << std::endl;
	PrintAugmentedMatrix(resultMatrix, finalVector, size); //TODO Working here

	//Stop console from closing
	int input;
	std::cin >> input;

	return 0;
}
//Andrew Sheridan
//Math 5610 
//Written in C++
#include "../../HW4/BackSubstitution/BackSubstitution/Matrix.h"
#include <iostream>

int main(void) {
	const unsigned size = 4;

	//Problem 1
	/*double** matrix = CreateDiagonallyDominantMatrix(size);
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
	PrintMatrix(L, size);*/


	//Problem 1
	double** matrix;
	double* vector;
	double* onesVector = CreateOnesVector(size);
	double** matrixCopy;
	double** resultMatrix;
	double* resultVector;

	for (int i = 10; i <= 160;  i *= 2) {
		matrix = CreateDiagonallyDominantMatrix(i);
		vector = CreateOnesVector(i);
		vector = VectorMatrixMultiply(matrix, vector, i);
		matrixCopy = CopyMatrix(matrix, i);
		ScaledLUFactorization(matrixCopy, vector, i);
		resultVector = BackSubstitution(matrixCopy, vector, i);

		/*std::cout << "Problem 5: Testing" << std::endl;
		std::cout << "Start matrix" << std::endl;
		PrintMatrix(matrix, i);
		std::cout << "Test vector " << std::endl;
		PrintVector(vector, i);
		std::cout << "Result of LU Factorization and Back Substitution " << std::endl;
		PrintAugmentedMatrix(matrixCopy, resultVector, i);*/

		std::cout << "Result of LU Factorization and Back Substitution with size " << i <<  std::endl;
		PrintVector(resultVector, i);
	}

	//Problem 2
	for (int i = 10; i <= 160; i *= 2) {
		matrix = CreateDiagonallyDominantMatrix(i);
		vector = CreateOnesVector(i);
		vector = VectorMatrixMultiply(matrix, vector, i);
		matrixCopy = CopyMatrix(matrix, i);
		matrixCopy = GaussianEliminationWithScaledPivoting(matrixCopy, vector, i);
		resultVector = BackSubstitution(matrixCopy, vector, i);

		/*std::cout << "Problem 5: Testing" << std::endl;
		std::cout << "Start matrix" << std::endl;
		PrintMatrix(matrix, i);
		std::cout << "Test vector " << std::endl;
		PrintVector(vector, i);
		std::cout << "Result of LU Factorization and Back Substitution " << std::endl;
		PrintAugmentedMatrix(matrixCopy, resultVector, i);*/

		std::cout << "Result of GE with Scaled Pivoting and Back Substitution with size " << i << std::endl;
		PrintVector(resultVector, i);
	}

	//Problem 3
	for (int i = 10; i <= 160; i *= 2) {
		matrix = CreateSymmetricMatrix(i);
		vector = CreateOnesVector(i);
		vector = VectorMatrixMultiply(matrix, vector, i);
		matrixCopy = CopyMatrix(matrix, i);
		GaussianElimination(matrixCopy, vector, i);
		resultVector = BackSubstitution(matrixCopy, vector, i);

		/*std::cout << "Problem 5: Testing" << std::endl;
		std::cout << "Start matrix" << std::endl;
		PrintMatrix(matrix, i);
		std::cout << "Test vector " << std::endl;
		PrintVector(vector, i);
		std::cout << "Result of LU Factorization and Back Substitution " << std::endl;
		PrintAugmentedMatrix(matrixCopy, resultVector, i);*/

		std::cout << "Result of Gaussian Elimination and Back Substitution with size " << i << std::endl;
		PrintVector(resultVector, i);
	}

	int input;
	std::cin >> input;

	return 0;
}

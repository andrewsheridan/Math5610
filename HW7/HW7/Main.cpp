#include "Matrix.h"
#include "Vector.h"
#include "MatrixOperations.h"
#include "MatrixFactory.h"
#include <iostream>

int main() {
	const unsigned size = 3;
	MatrixFactory MF;

	/*
	//Problem 5
	Matrix matrix(size);
	matrix.InitializeDiagonallyDominant();

	Vector vector(size);
	vector[0] = 1;
	vector[1] = 2;
	vector[2] = 3;

	Vector fivePartOne = matrix * vector;
	std::cout << "Problem 5: Testing" << std::endl;
	std::cout << "--------------5.1--------------" << std::endl;
	matrix.Print();
	vector.Print();
	fivePartOne.Print();

	std::cout << "--------------5.3--------------" << std::endl;
	Vector onesVector(size);
	onesVector.InitializeAllOnes();
	Vector fivePartThree = matrix * onesVector;

	onesVector.Print();
	fivePartThree.Print();

	std::cout << "--------------5.4--------------" << std::endl;
	for (int i = 2; i <= 4; i++) {
		Matrix matrix(i);
		matrix.InitializeDiagonallyDominant();
		Vector vector(i);
		vector.InitializeAllOnes();
		Vector product = matrix * vector;

		std::cout << "Starting System" << std::endl;
		matrix.PrintAugmented(product);
		Matrix matrixCopy = matrix;
		std::cout << "After Gaussian Elimination" << std::endl;
		GaussianElimination(matrixCopy, product);
		matrixCopy.PrintAugmented(product);
		Vector resultVector = BackSubstitution(matrixCopy, product);
		std::cout << "Result of GE and Back Substitution " << std::endl;
		matrixCopy.PrintAugmented(resultVector);
	}
	*/
	//Problem 1: Cholesky Decomposition
	int size1 = 4;
	Matrix matrix1 = MF.SPD(size1);
	Matrix cholesky1 = CholeskyDecomposition(matrix1);
	if (cholesky1 != NULL) {
		std::cout << "Our diagonally dominant test matrix and the result of computing the Cholesky Decomposition of the matrix." << std::endl;
		matrix1.Print();
		cholesky1.Print();

		Matrix transpose = cholesky1.Transpose();
		std::cout << "L Transpose" << std::endl;
		transpose.Print();
		std::cout << "Multiplying the Cholesky Decomposition by its transpose." << std::endl;
		Matrix result = cholesky1 * transpose;
		result.Print();
	}

	//Problem 2: 
	// Testing Positive Definiteness and Symmetry with Cholesky
	int size2 = 4;

	//Test with a diagonally dominant symmetric matrix. 
	Matrix matrix2a = MF.SPD(size2);
	Matrix cholesky2a = CholeskyDecomposition(matrix2a);
	std::cout << "Diagonally Dominant Symmetric Matrix Test" << std::endl;
	if (cholesky2a == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl;

	//Test with a diagonally dominant asymmetric matrix
	Matrix matrix2b = MF.DiagonallyDominant(size2, size2);
	Matrix cholesky2b = CholeskyDecomposition(matrix2b);
	std::cout << "Diagonally Dominant Matrix Test" << std::endl;
	if (cholesky2b == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl;

	//Test with a symmetric matrix
	Matrix matrix2c = MF.Symmetric(size2);
	Matrix cholesky2c = CholeskyDecomposition(matrix2c);
	std::cout << "Symmetric Matrix Test" << std::endl;
	if (cholesky2c == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl;

	//Test with a matrix
	Matrix matrix2d = MF.Random(size2, size2);
	Matrix cholesky2d = CholeskyDecomposition(matrix2d);
	std::cout << "Matrix Test" << std::endl;
	if (cholesky2d == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl << std::endl;
	
	int input;
	std::cin >> input;

	return 0;
}
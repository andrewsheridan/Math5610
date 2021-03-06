//Andrew Sheridan
//Math 5610 
//Written in C++

//main.cpp
#include "Matrix.h";
#include "Vector.h"
#include <iostream>

int main() 
{
	//Problem 1: Cholesky Decomposition
	int size1 = 4;
	double** matrix1 = CreateDiagonallyDominantSymmetricMatrix(size1);
	double** cholesky1 = CholeskyDecomposition(matrix1, size1);
	if (cholesky1 != NULL) {
		std::cout << "Our diagonally dominant test matrix and the result of computing the Cholesky Decomposition of the matrix." << std::endl;
		PrintMatrix(matrix1, size1);
		PrintMatrix(cholesky1, size1);

		double** transpose = Transpose(cholesky1, size1);
		std::cout << "L Transpose" << std::endl;
		PrintMatrix(transpose, size1);
		std::cout << "Multiplying the Cholesky Decomposition by its transpose." << std::endl;
		double** result = DotProduct(cholesky1, transpose, size1, size1, size1);
		PrintMatrix(result, size1);
	}

	//Problem 2: 
	// Testing Positive Definiteness and Symmetry with Cholesky
	int size2 = 4;

	//Test with a diagonally dominant symmetric matrix. 
	double** matrix2a = CreateDiagonallyDominantSymmetricMatrix(size2);
	double** cholesky2a = CholeskyDecomposition(matrix2a, size2);
	std::cout << "Diagonally Dominant Symmetric Matrix Test" << std::endl;
	if (cholesky2a == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl;

	//Test with a diagonally dominant asymmetric matrix
	double** matrix2b = CreateDiagonallyDominantMatrix(size2);
	double** cholesky2b = CholeskyDecomposition(matrix2b, size2);
	std::cout << "Diagonally Dominant Matrix Test" << std::endl;
	if (cholesky2b == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl;
	
	//Test with a symmetric matrix
	double** matrix2c = CreateSymmetricMatrix(size2);
	double** cholesky2c = CholeskyDecomposition(matrix2c, size2);
	std::cout << "Symmetric Matrix Test" << std::endl;
	if (cholesky2c == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl;

	//Test with a matrix
	double** matrix2d = CreateMatrix(size2);
	double** cholesky2d = CholeskyDecomposition(matrix2d, size2);
	std::cout << "Matrix Test" << std::endl;
	if (cholesky2d == NULL)
		std::cout << "This matrix is not symmetric and positive definite." << std::endl;
	else
		std::cout << "This matrix is symmetric and positive definite. " << std::endl << std::endl;

	//Problem 3: 1-Norm
	int size3 = 4;
	double** matrix3 = CreateMatrix(size3);
	matrix3[size3 - 2][size3 - 1] *= 10;
	double result3 = OneNorm(matrix3, size3);

	std::cout << "Our test matrix and the result of computing the 1-Norm of the matrix." << std::endl;
	PrintMatrix(matrix3, size3);
	std::cout << result3 << std::endl << std::endl;

	//Problem 4: Infinity-Norm
	int size4 = 4;
	double** matrix4 = CreateMatrix(size4);
	matrix4[size4 - 2][size4 - 1] *= 10;
	double result4 = InfinityNorm(matrix4, size4);

	std::cout << "Our test matrix and the result of computing the Infinity-Norm of the matrix." << std::endl;
	PrintMatrix(matrix4, size4);
	std::cout << result4 << std::endl << std::endl;


	//Problem 5: Condition Number
	int size5 = 5;
	double** matrix5 = CreateMatrixWithRange(size5, 10000, -10000);
	double result5 = ConditionNumber(matrix5, size5);

	std::cout << "Our test matrix and the result of computing condition number." << std::endl;
	PrintMatrix(matrix5, size5);
	std::cout << "Condition Number: " << result5 << std::endl << std::endl;

	//Problem 6: Condition Number Operation Counts
	for (int size6 = 5; size6 <= 160; size6 *= 2) {
		double** matrix6 = CreateMatrix(size6);
		double conditionNumber = ConditionNumber(matrix6, size6);
		//std::cout << "Condition Number: " << conditionNumber << std::endl << std::endl;
	}


	//Problem 7: Tridiagonal Matrices
	int size7 = 5;

	std::cout << "TESTING MINIFIED TRIDIAGONALS" << std::endl;
	std::cout << "Tridiagonal matrix before gaussian elimination" << std::endl;
	double** matrix7 = CreateMinifiedTridiagonal(size7);
	double** copy7 = CopyMatrix(matrix7, size7);
	double* vector7 = CreateOnesVector(size7);
	PrintMatrix(matrix7, size7, 3);
	PrintVector(vector7, size7);

	std::cout << "Tridiagonal matrix after gaussian elimination" << std::endl;
	TridiagonalElimination(copy7, vector7, size7);
	PrintMatrix(copy7, size7, 3);
	PrintVector(vector7, size7);

	std::cout << "Solution after back substitution" << std::endl;
	double* result7 = TridiagonalBackSubstitution(copy7, vector7, size7);
	PrintVector(result7, size7);

	std::cout << "COMPARING TO FULL TRIDIAGONALS MATRIX OPERATIONS" << std::endl;
	std::cout << "Tridiagonal matrix before gaussian elimination" << std::endl;
	double** matrix7b = CreateTridiagonalMatrix(size7);
	PrintMatrix(matrix7b, size7);
	double* vector7b = CreateOnesVector(size7);

	std::cout << "Tridiagonal matrix after gaussian elimination" << std::endl;
	GaussianElimination(matrix7b, vector7b, size7);
	PrintMatrix(matrix7b, size7);
	PrintVector(vector7b, size7);
	double* result7b = BackSubstitution(matrix7b, vector7b, size7);

	std::cout << "Solution without minification." << std::endl;
	PrintVector(result7b, size7);

	int input; 
	std::cin >> input;

	return 0;
}
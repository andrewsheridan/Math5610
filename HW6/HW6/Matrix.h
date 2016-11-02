//Andrew Sheridan
//Math 5610 
//Written in C++

//Matrix.h
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include "Vector.h"

#pragma region Initalization Operations
///Creates the identity matrix of size n
double** CreateIdentityMatrix(unsigned n) {
	double** matrix = new double*[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (unsigned j = 0; j < n; j++) {
			matrix[i][j] = 0;
		}
		matrix[i][i] = 1;
	}
	return matrix;
}


/// Generates a random  square matrix of size n.
// n: The size of the matrix
double** CreateMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (unsigned j = 0; j < n; j++) {
			matrix[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}

	return matrix;
}

/// Generates a random upper triangular square matrix of size n
// n: The size of the matrix
double** CreateUpperTriangularMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (unsigned j = 0; j < n; j++) {
			if (j < i) {
				matrix[i][j] = 0;
			}
			else {
				matrix[i][j] = dis(generator); //Assign entry in matrix to random number between 0 and 1
			}
		}
	}

	for (unsigned k = 0; k < n; k++) {
		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
	}

	return matrix;
}

/// Generates a random lower triangular square matrix of size n
// n: The size of the matrix
double** CreateLowerTriangularMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (unsigned j = 0; j < n; j++) {
			if (i < j) {
				matrix[i][j] = 0;
			}
			else {
				matrix[i][j] = dis(generator); //Assign entry in matrix to random number between 0 and 1
			}
		}
	}

	for (unsigned k = 0; k < n; k++) {
		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
	}

	return matrix;
}

/// Generates a random diagonally dominant square matrix of size n.
// n: The size of the matrix
double** CreateDiagonallyDominantMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n];
		for (unsigned j = 0; j < n; j++) {
			matrix[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}

	for (unsigned k = 0; k < n; k++) {
		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
	}

	return matrix;
}

/// Generates a random diagonally dominant square matrix of size n.
// n: The size of the matrix
double** CreateDiagonallyDominantSymmetricMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n]; //Must do this before our second loop, so that all rows are initialized.
	}
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = i; j < n; j++) {
			double value = dis(generator); //Assign each entry in matrix to random number between 0 and 1
			matrix[i][j] = value;
			matrix[j][i] = value;
		}
	}

	for (unsigned k = 0; k < n; k++) {
		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
	}

	return matrix;
}


/// Generates a symmetric square matrix of size n.
// n: The size of the matrix
double** CreateSymmetricMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[n]; //Must do this before our second loop, so that all rows are initialized.
	}
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = i; j < n; j++) {
			double value = dis(generator); //Assign each entry in matrix to random number between 0 and 1
			matrix[i][j] = value;
			matrix[j][i] = value;
		}
	}

	return matrix;
}

double** CreateMinifiedTridiagonalMatrix(unsigned n) {
	double** newMatrix = new double*[3];
	newMatrix[0] = new double[n - 1];
	newMatrix[1] = new double[n];
	newMatrix[2] = new double[n - 1];

	for (int i = 0; i < n - 1; i++) {
		newMatrix[0][i] = 1;
		newMatrix[1][i] = -2;
		newMatrix[2][i] = 1;
	}
	newMatrix[1][n - 1] = -2;
	
	return newMatrix;
}

double** CreateTridiagonalMatrix(unsigned n) {
	double** newMatrix = new double*[n];

	for (int i = 0; i < n - 1; i++) {
		newMatrix[i] = new double[n];
		newMatrix[i][i] = -2;
		newMatrix[i][i + 1] = 1;
		newMatrix[i + 1][i] = 1;
	}
	newMatrix[n - 1][n - 1] = -2;

	return newMatrix;
}
#pragma endregion

#pragma region Comparison Operations
///Returns a copy of the input nxn matrix
double** CopyMatrix(double** matrix, unsigned size) {
	double** result = new double*[size];
	for (unsigned i = 0; i < size; i++) {
		result[i] = new double[size];
		for (unsigned j = 0; j < size; j++) {
			result[i][j] = matrix[i][j];
		}
	}
	return result;
}
///Compares two nxn matrices, A and B. Returns true if they are identital, and false if they differ.
bool CompareMatrices(double** A, double** B, unsigned n) {
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			if (A[i][j] != B[i][j]) {
				return false;
			}
		}
	}
	return true;
}

//Checks to see if nxn matrix A is symmetric
bool IsMatrixSymmetric(double** A, unsigned n) {
	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j <= i; j++) {
			if (A[i][j] != A[j][i])
				return false;
		}
	}
	return true;
}


#pragma endregion

#pragma region HW4
/// Solves an nxn set of linear equations using back substitution
// A: The nxn upper-triangular coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double* BackSubstitution(double **A, double *b, unsigned n) {

	double* x;
	x = new double[n];
	try {
		x[n - 1] = b[n - 1] / A[n - 1][n - 1];
		for (int k = n - 2; k >= 0; k--) {
			x[k] = b[k];
			for (int i = k + 1; i < n; i++) {
				x[k] -= A[k][i] * x[i];
			}
			x[k] /= A[k][k];
		}
	}
	catch (std::exception& e)
	{
		std::cout << "These matrices are not the correct size." << std::endl;
		return new double[n];
	}

	return x;
}

/// Solves a set of linear equations using forward substitution
//a: the nxn lower-triangular coefficient matrix
//b: right-hand-side
//n: the size of the matrices
double* ForwardSubstitution(double** A, double* b, unsigned n) {
	double* x;
	x = new double[n];
	try {
		x[0] = b[0];
		for (int k = 0; k < n; k++) {
			x[k] = b[k];
			for (int j = 0; j < k; j++) {
				x[k] = x[k] - A[k][j] * x[j];
			}
			x[k] = x[k] / A[k][k];
		}
	}
	catch (std::exception& e)
	{
		std::cout << "These matrices are not the correct size." << std::endl;
		return new double[n];
	}

	return x;
}

/// Reduces an nxn matrix A and right-hand-side b to upper-triangular form using Gaussian Elimination
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double** GaussianElimination(double** A, double* b, unsigned n) {
	try {
		for (int k = 0; k < n; k++) {
			for (int i = k + 1; i < n; i++) {
				double factor = A[i][k] / A[k][k];
				for (int j = 0; j < n; j++) {
					A[i][j] = A[i][j] - factor*A[k][j];
				}
				b[i] = b[i] - factor*b[k];
			}
		}
	}
	catch (std::exception& e)
	{
		std::cout << "These matrices are not the correct size." << std::endl;
	}

	return A;
}
#pragma endregion

#pragma region HW5
/// Reduces an nxn matrix A and right-hand-side b to upper-triangular form using Gaussian Elimination
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double** GaussianEliminationWithScaledPivoting(double** A, double* b, unsigned n) {
	try {
		for (int k = 0; k < n; k++) {
			double* ratios = new double[n - k]; // New vector of size n - k to store the ratios
			for (int i = k; i < n; i++) {
				double rowMax = FindArrayMax(A[i], k, n);
				ratios[i - k] = rowMax / A[i][k];
			}
			int newPivot = FindMaxIndex(ratios, n - k) + k; //Find the best row for this iteration

			double* temp = A[k]; //
			A[k] = A[newPivot];  // Switch the current row with the best row for this iteration
			A[newPivot] = temp;  //

			double tempEntry = b[k];
			b[k] = b[newPivot];
			b[newPivot] = tempEntry;

			for (int i = k + 1; i < n; i++) {
				double factor = A[i][k] / A[k][k];
				for (int j = 0; j < n; j++) {
					A[i][j] = A[i][j] - factor*A[k][j];
				}
				b[i] = b[i] - factor*b[k];
			}
		}
	}
	catch (std::exception& e)
	{
		std::cout << "These matrices are not the correct size." << std::endl;
	}

	return A;
}


/// Finds the LU factorization of matrix A. A becomes the upper triangular matrix U, and the lower triangular matrix L is returned. 
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double** LUFactorization(double** A, double* b, unsigned n) {
	double** L = CreateIdentityMatrix(n);
	try {
		for (int k = 0; k < n; k++) {
			for (int i = k + 1; i < n; i++) {
				double factor = A[i][k] / A[k][k];
				L[i][k] = factor;
				for (int j = 0; j < n; j++) {
					A[i][j] = A[i][j] - factor*A[k][j];
				}
				b[i] = b[i] - factor*b[k];
			}
		}
	}
	catch (std::exception& e)
	{
		std::cout << "These matrices are not the correct size." << std::endl;
	}

	return L;
}

/// Finds the LU factorization of matrix A. A becomes the upper triangular matrix U, and the lower triangular matrix L is returned. 
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double** ScaledLUFactorization(double** A, double* b, unsigned n) {
	double** L = CreateIdentityMatrix(n);

	for (int k = 0; k < n; k++) {
		double* ratios = new double[n - k]; // New vector of size n - k to store the ratios
		for (int i = k; i < n; i++) {
			double rowMax = FindArrayMax(A[i], k, n);
			ratios[i - k] = rowMax / A[i][k];
		}
		int newPivot = FindMaxIndex(ratios, n - k) + k; //Find the best row for this iteration

		double* temp = A[k]; //
		A[k] = A[newPivot];  // Switch the current row with the best row for this iteration
		A[newPivot] = temp;  //

		double tempEntry = b[k];
		b[k] = b[newPivot];
		b[newPivot] = tempEntry;
		/*if (newPivot != k) {
		std::cout << "Exchanged rows " << k << " and " << newPivot << std::endl;
		}*/

		for (int i = k + 1; i < n; i++) {
			double factor = A[i][k] / A[k][k];
			L[i][k] = factor;
			for (int j = k + 1; j < n; j++) {
				A[i][j] = A[i][j] - factor*A[k][j];
			}
			A[i][k] = 0;
			b[i] = b[i] - factor*b[k];
		}
	}

	return L;
}

///Multiplies an nxn matrix A by the vector x
double* VectorMatrixMultiply(double** A, double* x, unsigned n) {
	double* result = new double[n];
	try {
		for (unsigned i = 0; i < n; i++) {
			result[i] = 0;
			for (unsigned j = 0; j < n; j++) {
				result[i] += A[i][j] * x[j];
			}
		}
		return result;
	}
	catch (std::exception& e) {
		std::cout << "These matrices are not the correct size." << std::endl;
	}
}
#pragma endregion

#pragma region HW6

double** Transpose(double** A, unsigned int n) {
	double** matrix = new double*[n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[j][i] = A[i][j];
		}
	}
	return matrix;
}

double** DotProduct(double**A, double** B, unsigned int m, unsigned int n, unsigned int p) {
	double** matrix = new double*[m];
	double** bTranspose = Transpose(B, n);
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double[p];
		for (unsigned j = 0; j < n; j++) {
			matrix[i][j] = DotProduct(A[i], bTranspose[j], n);
		}
	}
	return matrix;
}

///Computes the Cholesky Decomposition of an n by n matrix A
/// Returns NULL if the matrix is not SPD
double** CholeskyDecomposition(double** A, unsigned int n) {
	if (IsMatrixSymmetric(A, n) == false)
		return NULL;

	double** L = new double*[n]; //Initialize the new matrix
	for (int i = 0; i < n; i++) {  
		L[i] = new double[n];
		for (int j = 0; j < n; j++) {
			L[i][j] = 0;
		}
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < (i + 1); j++) {
			double entry = 0;
			for (int k = 0; k < j; k++) {
				entry += L[i][k] * L[j][k];
			}
			double sqrtValue = A[i][i] - entry;
			if (sqrtValue < 0)
				return NULL;
			L[i][j] = i == j ? std::sqrt(sqrtValue) : (1.0 / L[j][j] * (A[i][j] - entry));
		}
	}

	return L;
}

///Computes the 1-norm of an n by n matrix A
double OneNorm(double** A, unsigned int n) {
	double columnMax = 0;
	for (unsigned int i = 0; i < n; i++) {
		double columnSum = 0;
		for (unsigned int j = 0; j < n; j++) {
			columnSum += std::abs(A[i][j]);
		}
		if (columnSum > columnMax)
			columnMax = columnSum;
	}
	return columnMax;
}

///Computes the inifinity norm of an n by n matrix A
double InfinityNorm(double** A, unsigned int n) {
	double rowMax = 0;
	for (unsigned int j = 0; j < n; j++) {
		double rowSum = 0;
		for (unsigned int i = 0; i < n; i++) {
			rowSum += std::abs(A[i][j]);
		}
		if (rowSum > rowMax)
			rowMax = rowSum;
	}
	return rowMax;
}
#pragma endregion

#pragma region Printing

///Outputs an nxn matrix to the console
void PrintMatrix(double** matrix, unsigned size) {
	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = 0; j < size; j++) {
			std::cout << std::setw(13) << std::left << matrix[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}



///Outputs an augmented coefficient matrix to the console
void PrintAugmentedMatrix(double** matrix, double* vector, unsigned size) {
	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = 0; j < size; j++) {
			std::cout << std::setw(13) << std::left << matrix[i][j];
		}
		std::cout << "| " << vector[i] << std::endl;
	}
	std::cout << std::endl;
}
#pragma endregion


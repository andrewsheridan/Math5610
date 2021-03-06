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


///Outputs an n x m matrix to the console
void PrintMatrix(double** matrix, unsigned int n, unsigned int m) {
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < m; j++) {
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

/// Generates a random  square matrix of size n.
// n: The size of the matrix
double** CreateMatrixWithRange(unsigned n, int maxValue, int minValue) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(minValue, maxValue); //Desired distribution

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
/// n: The size of the matrix
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

double** CreateTridiagonalMatrix(unsigned n) {
	double** newMatrix = new double*[n];
	for (int i = 0; i < n; i++) {
		newMatrix[i] = new double[n];
		for (int j = 0; j < n; j++) {
			newMatrix[i][j] = 0;
		}
	}
	for (int i = 0; i < n - 1; i++) {
		newMatrix[i][i] = -2;
		newMatrix[i][i + 1] = 1;
		newMatrix[i + 1][i] = 1;
	}
	newMatrix[n - 1][n - 1] = -2;

	return newMatrix;
}

double** CreateMinifiedTridiagonal(unsigned n) {
	double** newMatrix = new double*[n];
	newMatrix[0] = new double[3];
	newMatrix[0][0] = 0;
	newMatrix[0][1] = -2;
	newMatrix[0][2] = 1;
	for (int i = 1; i < n - 1; i++) {
		newMatrix[i] = new double[3];
		newMatrix[i][0] = 1;
		newMatrix[i][1] = -2;
		newMatrix[i][2] = 1;
	}
	newMatrix[n - 1] = new double[3];
	newMatrix[n - 1][0] = 1;
	newMatrix[n - 1][1] = -2;
	newMatrix[n - 1][2] = 0;
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

///Checks to see if nxn matrix A is symmetric
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
/// Returns the transpose of the n by n matrix A
double** Transpose(double** A, unsigned int m, unsigned int n) {
	double** matrix = new double*[m];
	for (int i = 0; i < m; i++) {
		matrix[i] = new double[n];
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			matrix[j][i] = A[i][j];
		}
	}
	return matrix;
}

/// Computes the dot product of matrices A (m x n) and B (n x p)
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

			// Conditional assignment. If the entry is diagonal, assign to the sqare root of the previous value. 
			// Otherwise, Do computation for a nondiagonal entry.
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

///Computes the inifinity norm of an n by n matrix A
double InfinityNormOperations(double** A, unsigned int n, long& counter) {
	double rowMax = 0;
	for (unsigned int j = 0; j < n; j++) {
		double rowSum = 0;
		for (unsigned int i = 0; i < n; i++) {
			counter++;
			rowSum += std::abs(A[i][j]);
		}
		if (rowSum > rowMax)
			rowMax = rowSum;
		counter++;
	}
	return rowMax;
}

//Computes the inverse of n by n matrix A
double** Inverse(double** A, unsigned int n) {
	double** matrix = CreateIdentityMatrix(n);
	double ratio, a;
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i != j) {
				ratio = A[j][i] / A[i][i];
				for (k = 0; k < n; k++) {
					A[j][k] -= ratio * A[i][k];
				}
				for (k = 0; k < n; k++) {
					matrix[j][k] -= ratio * matrix[i][k];
				}
			}
		}
	}
	for (i = 0; i < n; i++) {
		a = A[i][i];
		for (j = 0; j < n; j++) {
			matrix[i][j] /= a;
		}
	}
	return matrix;
}

//Computes the inverse of n by n matrix A, incrementing the passed in counter
double** InverseOperations(double** A, unsigned int n, long& counter) {
	double** matrix = CreateIdentityMatrix(n);
	double ratio, a;
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			counter++;
			if (i != j) {
				ratio = A[j][i] / A[i][i];
				counter++;
				for (k = 0; k < n; k++) {
					A[j][k] -= ratio * A[i][k];
					counter++;
				}
				for (k = 0; k < n; k++) {
					matrix[j][k] -= ratio * matrix[i][k];
					counter++;
				}
			}
		}
	}
	for (i = 0; i < n; i++) {
		a = A[i][i];
		for (j = 0; j < n; j++) {
			counter++;
			matrix[i][j] /= a;
		}
	}
	return matrix;
}

//Estimates the condition number of n by n matrix A
double ConditionNumber(double** A, unsigned int n) {
	long counter = 0;
	double** aCopy = CopyMatrix(A, n);
	double** inverse = InverseOperations(aCopy, n, counter);

	double aNorm = InfinityNormOperations(A, n, counter);
	double inverseNorm = InfinityNormOperations(inverse, n, counter);

	std::cout << "The operations required to estimate CN for size " << n << ": " << counter << std::endl;
	return aNorm * inverseNorm;
}

//double** TridiagonalSimplification(double** A, unsigned int n) {
//	double** L = new double*[n];
//	for (int k = 0; k < n; k++) {
//		L[k] = new double[3];
//	}
//	for (int k = 0; k < n; k++) {
//		L[k + 1][k] = A[k + 1][k] / A[k][k];
//		for(int )
//	}
//}

//Finds the LU decomposition of n by n matrix A, with lower bandwitdh p and upper bandwith q
double** BandedMatrixLUDecomposition(double** A, unsigned int n, unsigned int p, unsigned int q) {
	double** L = CreateIdentityMatrix(n);
	for (int k = 0; k < n - 1; k++) {
		for (int i = k + 1; i < k + p; i++) {
			L[i][k] = A[i][k] / A[k][k]; 
			for (int j = k + 1; j < k + q; j++) {
				A[i][j] = A[i][j] - L[i][k] * A[k][j];
			}
		}
	}
	return L;
}

///Does Gaussian Elimination on a minified tridiagonal matrix
void TridiagonalElimination(double** A, unsigned int n) {
	for (unsigned int i = 0; i < n - 1; i++) {
		double factor = A[i + 1][0] / A[i][1];
		A[i + 1][0] = 0;
		A[i + 1][1] -= (A[i][2] * factor);
	}
}

///Does Gaussian Elimination on a minified tridiagonal matrix and right-hand-side b
void TridiagonalElimination(double** A, double* b, unsigned int n) {
	for (unsigned int i = 0; i < n - 1; i++) {
		double factor = A[i + 1][0] / A[i][1];
		A[i + 1][0] = 0;
		A[i + 1][1] -= (A[i][2] * factor);
		b[i + 1] -= b[i] * factor;
	}
}

///Does Back Substitution on a minified tridiagonal matrix and right-hand-side b
double* TridiagonalBackSubstitution(double** A, double* b, unsigned int n) {
	double* x = new double[n];
	x[n - 1] = b[n - 1] / A[n - 1][1];
	for (int k = n-2; k >= 0; k--) {
		x[k] = b[k];
		x[k] -= A[k][2] * x[k + 1];
		x[k] /= A[k][1];
	}
	return x;
}
#pragma endregion

#pragma region HW7

//double* LeastSquaresNormalEquations(double** A, double* b, unsigned int m, unsigned int n) {
//	double** At = Transpose(A, m, n);
//	double** B = DotProduct(A, At, n, m, n);
//	double* y = VectorMatrixMultiply(At, b, )
//
//	double** G = CholeskyDecomposition(B, n);
//	double** Gt = Transpose(G, n, n);
//
//	double** GEliminated = GaussianEliminationWithScaledPivoting(G, b, )
//}

#pragma endregion


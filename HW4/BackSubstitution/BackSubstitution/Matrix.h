//Andrew Sheridan
//Math 5610 
//Written in C++

//Matrix.h
#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>

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

///Outputs a size n vector to the console
void PrintVector(double* vector, unsigned size) {
	for (unsigned i = 0; i < size; i++) {
		std::cout << std::setw(13) << std::left << vector[i];
	}
	std::cout << std::endl << std::endl;
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

//Finds the entry in V (a size n vector) with the largest magnitude, starting with entry "start".
double FindArrayMax(double* V, unsigned start, unsigned size) {
	double max = 0;
	for (int i = start; i < size; i++) {
		double value = std::abs(V[i]);
		if (value > max) max = value;
	}
	return max;
}

//Finds the index of the value with the largest magnitude in a vector V with size n
int FindMaxIndex(double* V, unsigned n) {
	double max = 0;
	int index = -1;
	for (int i = 0; i < n; i++) {
		double value = std::abs(V[i]);
		if (value > max) {
			max = value;
			index = i;
		}
	}
	return index;
}

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
		if (newPivot != k) {
			std::cout << "Exchanged rows " << k << " and " << newPivot << std::endl;
		}

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

/// Generates a random diagonally dominant square matrix of size n.
// n: The size of the matrix
double** CreateDiagonallyDominantMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution
	
	double** matrix;
	matrix = new double *[n];
	for (unsigned i = 0; i < n; i++) {
		matrix[i] = new double [n];
		for (unsigned j = 0; j < n; j++) {
			matrix[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}

	for (unsigned k = 0; k < n; k++) {
		matrix[k][k] += 10*n; //Add 10*n to all diagonal entries
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

	for (unsigned k = 0; k < n; k++) {
		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
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

/// Generates a random vector of size n
// n: The size of the vector
double* CreateVector(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double* vector;
	vector = new double [n];
	for (unsigned i = 0; i < n; i++) {
		vector[i] = dis(generator); //Assign each entry in vector to random number between 0 and 1
	}

	return vector;
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

///Multiplies an nxn matrix A by another nxn matrix B
double** MatrixMultiply(double** A, double** B, unsigned n) {
	double** result = new double*[n];
	for (unsigned i = 0; i < n; i++) {
		result[i] = new double[n];
	}
	try {
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = 0; j < n; j++) {
				result[i][j] = A[i][j] * B[j][i];
			}
		}
		return result;
	}
	catch (std::exception& e) {
		std::cout << "These matrices are not the correct size." << std::endl;
	}
}

///Creates a vector of all ones of size n
double* CreateOnesVector(unsigned n) {
	double* vector = new double[n];
	for (unsigned i = 0; i < n; i++) {
		vector[i] = 1;
	}
	return vector;
}

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

///Returns a copy of the input vector
double* CopyVector(double* vector, unsigned size) {
	double* result = new double[size];
	for (unsigned i = 0; i < size; i++) {
		result[i] = vector[i];
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
//Andrew Sheridan
//Math 5610 
//Written in C++

//Matrix.h
#pragma once
#include <iostream>
#include <cmath>
#include <random>

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
			for(int j = 0; j < k; j++){
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
	double** matrix;
	matrix = new double*[n];
	try {
		for (int i = 0; i < n; i++) {
			matrix[i] = new double[n];
			for (int j = 0; j < n; j++) {
				
			}
		}

	}
	catch (std::exception& e)
	{
		std::cout << "These matrices are not the correct size." << std::endl;
		return new double*[n];
	}

	return matrix;
}

/// Generates a random square matrix of size n
// n: The size of the matrix
double** CreateMatrix(unsigned n) {
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

///Creates a vector of all ones of size n
double* CreateOnesVector(unsigned n) {
	double* vector = new double[n];
	for (unsigned i = 0; i < n; i++) {
		vector[i] = 1;
	}
	return vector;
}

///Outputs an nxn matrix to the console
void PrintMatrix(double** matrix, unsigned size) {
	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = 0; j < size; j++) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

///Outputs a size n vector to the console
void PrintVector(double* vector, unsigned size) {
	for (unsigned i = 0; i < size; i++) {
		std::cout << vector[i] << " ";
	}
	std::cout << std::endl << std::endl;
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

//bool checkArraySize(double *A, unsigned size) {
//	if (sizeof((A)) / sizeof((A[0])) == size)
//		return true;
//	else
//		return false;
//}
//
//bool checkNByNMatrixSize(double **A, unsigned size) {
//	if (sizeof((A)) / sizeof((A[0])) != size)
//		return false;
//	
//	for (unsigned i = 0; i < size; i++) {
//
//	}
//}
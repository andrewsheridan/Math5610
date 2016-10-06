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
		for (unsigned k = n - 2; k >= 0; k--) {
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
		return new double[];
	}

	return x;
}

/// Solves a set of linear equations using forward substitution
// A: The nxn lower-triangular coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double* forwardSubstitution(double** A, double* b, unsigned n) {


}

/// Reduces an nxn matrix A and right-hand-side b to upper-triangular form using Gaussian Elimination
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
double** GaussianElimination(double** A, double* b, unsigned n) {

}

/// Generates a random square matrix of size n
// n: The size of the matrix
double** GenerateMatrix(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution
	
	double** matrix;
	matrix = new double *[n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double [n];
		for (int j = 0; j < n; j++) {
			matrix[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}

	for (int k = 0; k < n; k++) {
		matrix[k][k] += 10*n; //Add 10*n to all diagonal entries
	}

	return matrix;
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
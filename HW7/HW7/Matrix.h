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

class Matrix{
public: 
	
	Matrix() = default;
	Matrix(unsigned size);
	Matrix(unsigned rowCount, unsigned columnCount);
	Matrix(const Matrix &m);
	//Matrix(Matrix& m);
	Matrix operator = (const Matrix& m);
	~Matrix();
	
	void InitializeIdentityMatrix();
	void InitializeRandom();
	void InitializeRange(double minValue, double maxValue);
	void InitializeDiagonallyDominant();
	
	void Print();
	void PrintAugmented(Vector v);

	Vector &operator[] (unsigned row) { return entries[row]; }
	friend bool operator == (const Matrix& A, const Matrix& B);
	friend bool operator != (const Matrix& A, const Matrix& B);
	friend Vector operator * (const Matrix& A, Vector& x);
	friend Matrix operator * (Matrix A, Matrix B);

	bool IsSymmetric();
	Matrix Transpose();
	double OneNorm();
	double InfinityNorm();

	unsigned GetRows() { return rows; }
	unsigned GetColumns() { return columns; }
	void SetRows(unsigned r) { rows = r; }
	void SetColumns(unsigned c) { columns = c; }

private:
	Vector* entries;

	unsigned rows;
	unsigned columns;

};

#pragma region Initalization Operations

///// Generates a random upper triangular square matrix of size n
//// n: The size of the matrix
//Matrix CreateUpperTriangularMatrix(unsigned n) {
//	std::mt19937 generator(123); //Random number generator
//	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution
//
//	Matrix matrix;
//	matrix = new double *[n];
//	for (unsigned i = 0; i < n; i++) {
//		matrix[i] = new double[n];
//		for (unsigned j = 0; j < n; j++) {
//			if (j < i) {
//				matrix[i][j] = 0;
//			}
//			else {
//				matrix[i][j] = dis(generator); //Assign entry in matrix to random number between 0 and 1
//			}
//		}
//	}
//
//	for (unsigned k = 0; k < n; k++) {
//		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
//	}
//
//	return matrix;
//}
//
///// Generates a random lower triangular square matrix of size n
//// n: The size of the matrix
//Matrix CreateLowerTriangularMatrix(unsigned n) {
//	std::mt19937 generator(123); //Random number generator
//	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution
//
//	Matrix matrix;
//	matrix = new double *[n];
//	for (unsigned i = 0; i < n; i++) {
//		matrix[i] = new double[n];
//		for (unsigned j = 0; j < n; j++) {
//			if (i < j) {
//				matrix[i][j] = 0;
//			}
//			else {
//				matrix[i][j] = dis(generator); //Assign entry in matrix to random number between 0 and 1
//			}
//		}
//	}
//
//	for (unsigned k = 0; k < n; k++) {
//		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
//	}
//
//	return matrix;
//}
//
///// Generates a random diagonally dominant square matrix of size n.
///// n: The size of the matrix
//Matrix CreateDiagonallyDominantSymmetricMatrix(unsigned n) {
//	std::mt19937 generator(123); //Random number generator
//	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution
//
//	Matrix matrix;
//	matrix = new double *[n];
//	for (unsigned i = 0; i < n; i++) {
//		matrix[i] = new double[n]; //Must do this before our second loop, so that all rows are initialized.
//	}
//	for (unsigned i = 0; i < n; i++) {
//		for (unsigned j = i; j < n; j++) {
//			double value = dis(generator); //Assign each entry in matrix to random number between 0 and 1
//			matrix[i][j] = value;
//			matrix[j][i] = value;
//		}
//	}
//
//	for (unsigned k = 0; k < n; k++) {
//		matrix[k][k] += 10 * n; //Add 10*n to all diagonal entries
//	}
//
//	return matrix;
//}
//
//
///// Generates a symmetric square matrix of size n.
//// n: The size of the matrix
//Matrix CreateSymmetricMatrix(unsigned n) {
//	std::mt19937 generator(123); //Random number generator
//	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution
//
//	Matrix matrix;
//	matrix = new double *[n];
//	for (unsigned i = 0; i < n; i++) {
//		matrix[i] = new double[n]; //Must do this before our second loop, so that all rows are initialized.
//	}
//	for (unsigned i = 0; i < n; i++) {
//		for (unsigned j = i; j < n; j++) {
//			double value = dis(generator); //Assign each entry in matrix to random number between 0 and 1
//			matrix[i][j] = value;
//			matrix[j][i] = value;
//		}
//	}
//
//	return matrix;
//}
//
//Matrix CreateTridiagonalMatrix(unsigned n) {
//	Matrix newMatrix = new Vector[n];
//	for (unsigned i = 0; i < n; i++) {
//		newMatrix[i] = new double[n];
//		for (unsigned j = 0; j < n; j++) {
//			newMatrix[i][j] = 0;
//		}
//	}
//	for (unsigned i = 0; i < n - 1; i++) {
//		newMatrix[i][i] = -2;
//		newMatrix[i][i + 1] = 1;
//		newMatrix[i + 1][i] = 1;
//	}
//	newMatrix[n - 1][n - 1] = -2;
//
//	return newMatrix;
//}
//
//Matrix CreateMinifiedTridiagonal(unsigned n) {
//	Matrix newMatrix = new Vector[n];
//	newMatrix[0] = new double[3];
//	newMatrix[0][0] = 0;
//	newMatrix[0][1] = -2;
//	newMatrix[0][2] = 1;
//	for (unsigned i = 1; i < n - 1; i++) {
//		newMatrix[i] = new double[3];
//		newMatrix[i][0] = 1;
//		newMatrix[i][1] = -2;
//		newMatrix[i][2] = 1;
//	}
//	newMatrix[n - 1] = new double[3];
//	newMatrix[n - 1][0] = 1;
//	newMatrix[n - 1][1] = -2;
//	newMatrix[n - 1][2] = 0;
//	return newMatrix;
//}
//#pragma endregion


//
/////Computes the inifinity norm of an n by n matrix A
//double InfinityNormOperations(Matrix A, unsigned n, long& counter) {
//	double rowMax = 0;
//	for (unsigned j = 0; j < n; j++) {
//		double rowSum = 0;
//		for (unsigned i = 0; i < n; i++) {
//			counter++;
//			rowSum += std::abs(A[i][j]);
//		}
//		if (rowSum > rowMax)
//			rowMax = rowSum;
//		counter++;
//	}
//	return rowMax;
//}
//
//
////Computes the inverse of n by n matrix A, incrementing the passed in counter
//Matrix InverseOperations(Matrix A, unsigned n, long& counter) {
//	Matrix matrix = CreateIdentityMatrix(n);
//	double ratio, a;
//	unsigned i, j, k;
//	for (i = 0; i < n; i++) {
//		for (j = 0; j < n; j++) {
//			counter++;
//			if (i != j) {
//				ratio = A[j][i] / A[i][i];
//				counter++;
//				for (k = 0; k < n; k++) {
//					A[j][k] -= ratio * A[i][k];
//					counter++;
//				}
//				for (k = 0; k < n; k++) {
//					matrix[j][k] -= ratio * matrix[i][k];
//					counter++;
//				}
//			}
//		}
//	}
//	for (i = 0; i < n; i++) {
//		a = A[i][i];
//		for (j = 0; j < n; j++) {
//			counter++;
//			matrix[i][j] /= a;
//		}
//	}
//	return matrix;
//}
//
////Estimates the condition number of n by n matrix A
//double ConditionNumber(Matrix A, unsigned n) {
//	long counter = 0;
//	Matrix aCopy = CopyMatrix(A, n);
//	Matrix inverse = InverseOperations(aCopy, n, counter);
//
//	double aNorm = InfinityNormOperations(A, n, counter);
//	double inverseNorm = InfinityNormOperations(inverse, n, counter);
//
//	std::cout << "The operations required to estimate CN for size " << n << ": " << counter << std::endl;
//	return aNorm * inverseNorm;
//}
//
////Matrix TridiagonalSimplification(Matrix A, unsigned n) {
////	Matrix L = new Vector[n];
////	for (unsigned k = 0; k < n; k++) {
////		L[k] = new double[3];
////	}
////	for (unsigned k = 0; k < n; k++) {
////		L[k + 1][k] = A[k + 1][k] / A[k][k];
////		for(unsigned )
////	}
////}
//
////Finds the LU decomposition of n by n matrix A, with lower bandwitdh p and upper bandwith q
//Matrix BandedMatrixLUDecomposition(Matrix A, unsigned n, unsigned p, unsigned q) {
//	Matrix L = CreateIdentityMatrix(n);
//	for (unsigned k = 0; k < n - 1; k++) {
//		for (unsigned i = k + 1; i < k + p; i++) {
//			L[i][k] = A[i][k] / A[k][k];
//			for (unsigned j = k + 1; j < k + q; j++) {
//				A[i][j] = A[i][j] - L[i][k] * A[k][j];
//			}
//		}
//	}
//	return L;
//}
//
/////Does Gaussian Elimination on a minified tridiagonal matrix
//void TridiagonalElimination(Matrix A, unsigned n) {
//	for (unsigned i = 0; i < n - 1; i++) {
//		double factor = A[i + 1][0] / A[i][1];
//		A[i + 1][0] = 0;
//		A[i + 1][1] -= (A[i][2] * factor);
//	}
//}
//
/////Does Gaussian Elimination on a minified tridiagonal matrix and right-hand-side b
//void TridiagonalElimination(Matrix A, Vector b, unsigned n) {
//	for (unsigned i = 0; i < n - 1; i++) {
//		double factor = A[i + 1][0] / A[i][1];
//		A[i + 1][0] = 0;
//		A[i + 1][1] -= (A[i][2] * factor);
//		b[i + 1] -= b[i] * factor;
//	}
//}
//
/////Does Back Substitution on a minified tridiagonal matrix and right-hand-side b
//Vector TridiagonalBackSubstitution(Matrix A, Vector b, unsigned n) {
//	Vector x = new double[n];
//	x[n - 1] = b[n - 1] / A[n - 1][1];
//	for (unsigned k = n - 2; k >= 0; k--) {
//		x[k] = b[k];
//		x[k] -= A[k][2] * x[k + 1];
//		x[k] /= A[k][1];
//	}
//	return x;
//}
//#pragma endregion
//

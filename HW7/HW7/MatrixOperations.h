#pragma once
#include "Matrix.h"
#include "Vector.h"

#pragma region extra functions
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

#pragma endregion

#pragma region HW4
/// Solves a set of linear equations using back substitution
/// Does not reduce matrix A
//A: The matrix to be reduced
// b: Right-Hand-Side
Vector BackSubstitution(Matrix A, Vector b) {
	if (A.rows != b.size) return NULL;

	Vector x(b.size);

	x[b.size - 1] = b[b.size - 1] / A[A.rows - 1][A.columns - 1];
	for (int k = A.rows - 2; k >= 0; k--) {
		x[k] = b[k];
		for (unsigned i = k + 1; i < A.columns; i++) {
			x[k] -= A[k][i] * x[i];
		}
		x[k] /= A[k][k];
	}

	return x;
}

/// Solves a set of linear equations using forward substitution
/// Does not reduce matrix A
//A: The matrix to be reduced
//b: right-hand-side
Vector ForwardSubstitution(Matrix A, Vector b) {
	if (A.rows != b.size) return NULL;

	Vector x(b.size);

	x[0] = b[0];
	for (unsigned k = 0; k < A.rows; k++) {
		x[k] = b[k];
		for (unsigned j = 0; j < A.columns; j++) {
			x[k] = x[k] - A[k][j] * x[j];
		}
		x[k] = x[k] / A[k][k];
	}

	return x;
}

/// Reduces a matrix right-hand-side b to upper-triangular form using Gaussian Elimination
//A: The matrix to be reduced
// b: Right-Hand-Side
void GaussianElimination(Matrix& A, Vector& b) {
	for (unsigned k = 0; k < A.rows; k++) {
		for (unsigned i = k + 1; i < A.rows; i++) {
			double factor = A[i][k] / A[k][k];
			for (unsigned j = 0; j < A.columns; j++) {
				A[i][j] = A[i][j] - factor*A[k][j];
			}
			A[i][k] = 0;
			b[i] = b[i] - factor*b[k];
		}
	}
}
#pragma endregion

#pragma region HW5
/// Reduces an matrix A and right-hand-side b to upper-triangular form using Gaussian Elimination
// A: The coefficient matrix
// b: Right-Hand-Side
void GaussianEliminationWithScaledPivoting(Matrix A, Vector b) {
	if (A.rows != b.size) return;

	int n = b.size;

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


/// Finds the LU factorization of matrix A. A becomes the upper triangular matrix U, and the lower triangular matrix L is returned. 
// A: The nxn coefficient matrix
// b: Right-Hand-Side
Matrix LUFactorization(Matrix A, Vector b) {
	if (A.rows != b.size) return NULL;

	Matrix L(A.rows, A.columns);
	L.InitializeIdentityMatrix();

	for (int k = 0; k < A.rows; k++) {
		for (int i = k + 1; i < A.rows; i++) {
			double factor = A[i][k] / A[k][k];
			L[i][k] = factor;
			for (int j = 0; j < A.columns; j++) {
				A[i][j] = A[i][j] - factor*A[k][j];
			}
			b[i] = b[i] - factor*b[k];
		}
	}

	return L;
}

/// Finds the LU factorization of matrix A. A becomes the upper triangular matrix U, and the lower triangular matrix L is returned. 
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
Matrix ScaledLUFactorization(Matrix A, Vector b) {
	if (A.rows != b.size) return NULL;
	Matrix L(A.rows, A.columns);
	L.InitializeIdentityMatrix();

	for (unsigned k = 0; k < A.rows; k++) {
		Vector ratios(A.rows - k); // New vector of size A.rows - k to store the ratios
		for (unsigned i = k; i < A.rows; i++) {
			double rowMax = FindArrayMax(A[i], k, A.columns);
			ratios[i - k] = rowMax / A[i][k];
		}
		unsigned newPivot = ratios.FindMaxIndex() + k; //Find the best row for this iteration

		double* temp = A[k]; //
		A[k] = A[newPivot];  // Switch the current row with the best row for this iteration
		A[newPivot] = temp;  //

		double tempEntry = b[k];
		b[k] = b[newPivot];
		b[newPivot] = tempEntry;
		/*if (newPivot != k) {
		std::cout << "Exchanged rows " << k << " and " << newPivot << std::endl;
		}*/

		for (unsigned i = k + 1; i < A.rows; i++) {
			double factor = A[i][k] / A[k][k];
			L[i][k] = factor;
			for (unsigned j = k + 1; j < A.columns; j++) {
				A[i][j] = A[i][j] - factor*A[k][j];
			}
			A[i][k] = 0;
			b[i] = b[i] - factor*b[k];
		}
	}

	return L;
}

#pragma endregion

///Computes the Cholesky Decomposition of an n by n matrix A
/// Returns NULL if the matrix is not SPD
Matrix CholeskyDecomposition(Matrix& A) {
	if (A.IsSymmetric() == false)
		return NULL;

	Matrix L(A.rows, A.columns); //Initialize the new matrix
	for (unsigned i = 0; i < A.rows; i++) {
		for (unsigned j = 0; j < A.columns; j++) {
			L[i][j] = 0;
		}
	}

	for (unsigned i = 0; i < A.rows; i++) {
		for (unsigned j = 0; j < (i + 1); j++) {
			double entry = 0;
			for (unsigned k = 0; k < j; k++) {
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

//Computes the inverse of matrix A
Matrix Inverse(Matrix A) {
	Matrix matrix(A.columns, A.rows);
	matrix.InitializeIdentityMatrix();
	double ratio, a;
	unsigned i, j, k;
	for (i = 0; i < A.rows; i++) {
		for (j = 0; j < A.columns; j++) {
			if (i != j) {
				ratio = A[j][i] / A[i][i];
				for (k = 0; k < A.rows; k++) {
					A[j][k] -= ratio * A[i][k];
				}
				for (k = 0; k < A.rows; k++) {
					matrix[j][k] -= ratio * matrix[i][k];
				}
			}
		}
	}
	for (i = 0; i < A.rows; i++) {
		a = A[i][i];
		for (j = 0; j < A.columns; j++) {
			matrix[i][j] /= a;
		}
	}
	return matrix;
}

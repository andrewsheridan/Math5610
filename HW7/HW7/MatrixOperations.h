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
	if (A.GetRows() != b.GetSize()) return NULL;

	Vector x(b.GetSize());

	for(int i = b.GetSize() - 1; i >= 0; i--)
	{
		x[i] = b[i];
		for (int j = i + 1; j < b.GetSize(); j++) {
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
		x.Print();
	}
	return x;
}

/// Solves a set of linear equations using forward substitution
/// Does not reduce matrix A
//A: The matrix to be reduced
//b: right-hand-side
Vector ForwardSubstitution(Matrix A, Vector b) {
	if (A.GetRows() != b.GetSize()) return NULL;

	Vector x(b.GetSize());

	x[0] = b[0];
	for (unsigned k = 0; k < A.GetRows(); k++) {
		x[k] = b[k];
		for (unsigned j = 0; j < A.GetColumns(); j++) {
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
	for (unsigned k = 0; k < A.GetRows(); k++) {
		for (unsigned i = k + 1; i < A.GetRows(); i++) {
			double factor = A[i][k] / A[k][k];
			for (unsigned j = 0; j < A.GetColumns(); j++) {
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
	if (A.GetRows() != b.GetSize()) return;

	int n = b.GetSize();

	for (int k = 0; k < n; k++) {
		double* ratios = new double[n - k]; // New vector of size n - k to store the ratios
		for (int i = k; i < n; i++) {
			//double rowMax = FindArrayMax(A[i], k, n);
			double rowMax = A[i].FindMaxMagnitudeStartingAt(k);
			ratios[i - k] = rowMax / A[i][k];
		}
		int newPivot = FindMaxIndex(ratios, n - k) + k; //Find the best row for this iteration

		Vector temp = A[k]; //
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
	if (A.GetRows() != b.GetSize()) return NULL;

	Matrix L(A.GetRows(), A.GetColumns());
	L.InitializeIdentityMatrix();

	for (int k = 0; k < A.GetRows(); k++) {
		for (int i = k + 1; i < A.GetRows(); i++) {
			double factor = A[i][k] / A[k][k];
			L[i][k] = factor;
			for (int j = 0; j < A.GetColumns(); j++) {
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
	if (A.GetRows() != b.GetSize()) return NULL;
	Matrix L(A.GetRows(), A.GetColumns());
	L.InitializeIdentityMatrix();

	for (unsigned k = 0; k < A.GetRows(); k++) {
		Vector ratios(A.GetRows() - k); // New vector of size A.GetRows() - k to store the ratios
		for (unsigned i = k; i < A.GetRows(); i++) {
			//double rowMax = FindArrayMax(A[i], k, A.GetColumns());
			double rowMax = A[i].FindMaxMagnitudeStartingAt(k);
			ratios[i - k] = rowMax / A[i][k];
		}
		unsigned newPivot = ratios.FindMaxIndex() + k; //Find the best row for this iteration

		Vector temp = A[k]; //
		A[k] = A[newPivot];  // Switch the current row with the best row for this iteration
		A[newPivot] = temp;  //

		double tempEntry = b[k];
		b[k] = b[newPivot];
		b[newPivot] = tempEntry;
		/*if (newPivot != k) {
		std::cout << "Exchanged rows " << k << " and " << newPivot << std::endl;
		}*/

		for (unsigned i = k + 1; i < A.GetRows(); i++) {
			double factor = A[i][k] / A[k][k];
			L[i][k] = factor;
			for (unsigned j = k + 1; j < A.GetColumns(); j++) {
				A[i][j] = A[i][j] - factor*A[k][j];
			}
			A[i][k] = 0;
			b[i] = b[i] - factor*b[k];
		}
	}

	return L;
}

#pragma endregion

#pragma region HW6
///Computes the Cholesky Decomposition of an n by n matrix A
/// Returns NULL if the matrix is not SPD
Matrix CholeskyDecomposition(Matrix& A) {
	if (A.IsSymmetric() == false)
		return NULL;

	Matrix L(A.GetRows(), A.GetColumns()); //Initialize the new matrix
	for (unsigned i = 0; i < A.GetRows(); i++) {
		for (unsigned j = 0; j < A.GetColumns(); j++) {
			L[i][j] = 0;
		}
	}

	for (unsigned i = 0; i < A.GetRows(); i++) {
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
	Matrix matrix(A.GetColumns(), A.GetRows());
	matrix.InitializeIdentityMatrix();
	double ratio, a;
	unsigned i, j, k;
	for (i = 0; i < A.GetRows(); i++) {
		for (j = 0; j < A.GetColumns(); j++) {
			if (i != j) {
				ratio = A[j][i] / A[i][i];
				for (k = 0; k < A.GetRows(); k++) {
					A[j][k] -= ratio * A[i][k];
				}
				for (k = 0; k < A.GetRows(); k++) {
					matrix[j][k] -= ratio * matrix[i][k];
				}
			}
		}
	}
	for (i = 0; i < A.GetRows(); i++) {
		a = A[i][i];
		for (j = 0; j < A.GetColumns(); j++) {
			matrix[i][j] /= a;
		}
	}
	return matrix;
}

#pragma endregion

#pragma region HW7

//Todo: Check this.
Vector LeastSquares(Matrix A, Vector b) {
	Matrix AT = A.Transpose();
	Matrix B = AT * A;
	Vector y = AT * b;
	std::cout << "A:" << std::endl;
	A.Print();
	AT.Print();
	std::cout << "B: " << std::endl;
	B.Print();
	B.Transpose().Print();
	y.Print();

	Vector X = BackSubstitution(B, y);
	std::cout << "X" << std::endl;
	X.Print();
	Matrix G = CholeskyDecomposition(B);
	std::cout << "Cholesky: " << std::endl;
	G.Print();

	Vector z = BackSubstitution(G, y);
	std::cout << "Z:" << std::endl;
	z.Print();
	Vector x = BackSubstitution(G.Transpose(), z);
	std::cout << "X: " << std::endl;
	x.Print();
	return x;
}

double* CopyArray(double* v, unsigned n) {
	double* newArray = new double[n];
	for (int i = 0; i < n; i++) {
		newArray[i] = v[i];
	}
	return newArray;
}

void MultiplyByConstant(double* v, unsigned n, double multiplier) {
	for (int i = 0; i < n; i++) {
		v[i] *= multiplier;
	}
}

//Subtracts every entry in a by the corresponding entry in b
double* VectorSubtraction(double*a, double* b, unsigned size) {
	double* returnValue = new double[size];
	for (int i = 0; i < size; i++) {
		returnValue[i] = a[i] - b[i];
	}
	return returnValue;
}

Vector GramSchmidt(Matrix A) {
	if (A.GetRows() != A.GetColumns()) return NULL;
	Matrix q(A.GetRows());
	Matrix r(A.GetRows());
	for (int j = 0; j < A.GetRows(); j++) {
		q[j] = A[j];
		for (int i = 0; i < j - 1; j++) {
			r[i][j] = q[j] * q[i];
			q[j] = q[j] - (q[i] * r[i][j]);
		}

		//TODO: write vector norms.
		//r[j][j] = norm
		// qj = qj/ rjj
	}
}
#pragma endregion
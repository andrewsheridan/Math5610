//Andrew Sheridan
//Math 5610 
//Written in C++
//MatrixOperations.h

#pragma once
#include "Matrix.h"
#include "Vector.h"
#include "MatrixFactory.h"
#include <iostream>
#include <fstream>
#include <cmath>

#pragma region Basic Vector Operations
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
/// Note: Does not reduce matrix A
//A: The upper-triangular matrix
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
	}
	return x;
}

/// Solves a set of linear equations using forward substitution
/// Does not reduce matrix A
//A: The lower-triangular matrix
//b: right-hand-side
Vector ForwardSubstitution(Matrix A, Vector b) {
	if (A.GetRows() != b.GetSize()) return NULL;

	Vector x(b.GetSize());

	x[0] = b[0];
	for (unsigned i = 0; i < A.GetRows(); i++) {
		x[i] = b[i];
		for (unsigned j = 0; j < i; j++) {
			x[i] = x[i] - (A[i][j] * x[j]);
		}
		x[i] = x[i] / A[i][i];
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
Matrix GaussianEliminationWithScaledPivoting(Matrix A, Vector& b) {
	if (A.GetRows() != b.GetSize()) return NULL;

	int n = b.GetSize();

	for (int k = 0; k < n; k++) {
		double* ratios = new double[n - k]; // New vector of size n - k to store the ratios
		for (int i = k; i < n; i++) {
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
			b[i] = b[i] - (factor*b[k]);
		}
	}
	return A;
}


/// Finds the LU factorization of matrix A.
/// RHS b will be modified
// A: The nxn coefficient matrix
// b: Right-Hand-Side
Matrix* LUFactorization(Matrix A, Vector& b) {
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

	Matrix* LU = new Matrix[2]{ L, A };
	return LU;
}

/// Finds the LU factorization of matrix A. A becomes the upper triangular matrix U, and the lower triangular matrix L is returned. 
// A: The nxn coefficient matrix
// b: Right-Hand-Side
// n: The size of the matrices
Matrix* ScaledLUFactorization(Matrix A, Vector b) {
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

	return new Matrix[2]{ L, A };
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

/// A Least Squares algorithm via Normal Equations
/// Requires a matrix A and a vector b
Vector LeastSquares(Matrix A, Vector b) {
	Matrix AT = A.Transpose();
	Matrix B = AT * A;
	Vector y = AT * b;

	std::cout << "B: " << std::endl;
	B.Print();
	std::cout << "Y: " << std::endl;
	y.Print();

	Matrix G = CholeskyDecomposition(B);
	Vector z = ForwardSubstitution(G, y);
	Vector x = BackSubstitution(G.Transpose(), z);

	return x;
}

///Computes the QR factorization of Matrix A
///Returns a pair of matrices in an array. The first is Q, the second, R. 
Matrix* GramSchmidt(Matrix A) {
	if (A.GetRows() != A.GetColumns()) return NULL;

	Matrix r(A.GetRows(), A.GetColumns());
	Matrix q(A.GetRows(), A.GetColumns());

	for (int k = 0; k < A.GetRows(); k++) {
		r[k][k] = 0; 
		for (int i = 0; i < A.GetRows(); i++)
			r[k][k] = r[k][k] + A[i][k] * A[i][k];

		r[k][k] = sqrt(r[k][k]);

		for (int i = 0; i < A.GetRows(); i++) 
			q[i][k] = A[i][k] / r[k][k];

		for (int j = k + 1; j < A.GetColumns(); j++) {
			r[k][j] = 0;
			for (int i = 0; i < A.GetRows(); i++) 
				r[k][j] += q[i][k] * A[i][j];

			for (int i = 0; i < A.GetRows(); i++)
				A[i][j] = A[i][j] - r[k][j] * q[i][k];
		}
	}
	return new Matrix[2]{ q, r };
}

///Creates a vector of data which can be used to test least squares methods
///Takes two vectors, t and b, and an int n (the size of the sample space)
Vector LSFit(Vector t, Vector b, int n) {
	int m = t.GetSize();

	Matrix A = MatrixFactory::Instance()->Ones(m, n);
	for (int j = 0; j < n - 1; j++) {
		for (int i = 0; i < m; i++) {
			A[i][j + 1] =  A[i][j] * t[i];
		}
	}
	A.Print();
	Matrix AT = A.Transpose();
	Matrix B = AT * A;
	Vector y = AT * b;

	return B / y;
}

#pragma endregion

#pragma region HW9
///Solves the system of equations using Jacobi Iteration
///A: The Matrix
///x0: The initial guess
///b: The Right-Hand-Side
Vector JacobiIteration(Matrix A, Vector x0, Vector b, int maxIterations, double tolerance) {
	int iterations = 0;
	int n = A.GetRows();
	Vector newX(x0);
	double error = 10 * tolerance;
	while(iterations < maxIterations && tolerance < error){
		Vector oldX(newX);
		for (int i = 0; i < n;  i++) {
			newX[i] = b[i];
			for (int j = 0; j < i; j++) {
				newX[i] = newX[i] - A[i][j] * oldX[j];
			}
			for (int j = i + 1; j < n; j++) {
				newX[i] = newX[i] - A[i][j] * oldX[j];
			}
			newX[i] = newX[i] / A[i][i];
			error = (oldX - newX).L2Norm();
			iterations++;
		}
	}
	return newX;
}

///Solves a matrix and RHS using Conjugate Gradient Method
///A : The matrix to be solved
///b : The RHS vector
///x0 : The initial guess vector
///tol : The tolerance of our method
Vector ConjugateGradient(Matrix A, Vector b, Vector x0, double tol) {
	Vector rk = b - (A * x0);
	double dk = rk * rk;
	double bd = b * b;
	int k = 0;
	Vector pk = rk;
	Vector xk = x0;
	while (dk > tol * tol * bd) {
		Vector sk = A * pk;
		double ak = dk / (pk * sk);
		Vector xkp1 = xk + (ak * pk);
		Vector rkp1 = rk - (ak * sk);
		double dkp1 = rkp1 * rkp1;
		Vector pkp1 = rkp1 + ((dkp1 / dk) * pk);
		k++;
		//All values have been computed, set all kth values to equal the k+1 value in preparation for next iteration
		xk = xkp1;
		rk = rkp1;
		pk = pkp1;
		dk = dkp1;
	}
	return xk;
}
#pragma endregion

#pragma region HW10 

Vector SolveSystem(Matrix A, Vector b) {
	GaussianElimination(A, b);
	return BackSubstitution(A, b);
}

///Finds an approximation of the largest Eigenvalue of a matrix
///A : The square matrix
///x0 : The initial guess vector
///tol : The tolerance of the algorithm
///maxIter : The maximum number of iterations to be executed by the method
double PowerMethod(Matrix A, Vector x0, double tol, int maxIter) {
	double error = 10 * tol;
	int k = 0; 
	Vector y = A * x0;
	Vector xk = x0;
	double lambda_k = 0;

	//Output for problem 10.2. Comment out if not needed.
	std::ofstream output("powerMethod.txt");
	output << "Iterations \t Error " << std::endl;
	
	while (error > tol && k < maxIter) {
		Vector xkp1 = y / y.L2Norm();
		y = A * xkp1;
		double lambda_kp1 = xkp1 * y;
		error = abs(lambda_kp1 - lambda_k);
		/*std::cout << "lambda_k: " << lambda_k << std::endl;
		std::cout << "lambda_k+1: " << lambda_kp1 << std::endl;*/

		//Output for problem 10.2. Comment out if not needed.
		output << k << "\t" << error << std::endl;
			
		lambda_k = lambda_kp1;
		k++;
	}

	output.close();

	/*std::cout << "Eigenvector approximation: " << std::endl;
	y.Print();*/

	return lambda_k;
}

///Finds an approximation of the smallest Eigenvalue of a matrix
///A : The square matrix
///x0 : The initial guess vector
///tol : The tolerance of the algorithm
///maxIter : The maximum number of iterations to be executed by the method
double InversePowerMethod(Matrix A, Vector x0, double tol, int maxIter) {

	//Output for problem 10.5. Comment out if not needed.
	std::ofstream output("inversePowerMethod.txt");
	output << "Iterations \t Error " << std::endl;

	double error = 10 * tol;
	int k = 0;
	Matrix* LU = LUFactorization(A, x0);
	Matrix U = LU[1];
	Matrix L = LU[0];
	Vector y = BackSubstitution(U, x0);// Solve for y by doing back substitution.
	double lambda_x = 0;
	while (error > tol && k < maxIter) {
		Vector x = y / y.L2Norm();
		y = SolveSystem(A, x); //Does Gaussian Elimination then Back Substitution
		double lambda_xp1 = x * y;
		error = abs(lambda_xp1 - lambda_x);
		/*std::cout << "lambda_k: " << lambda_x << std::endl;
		std::cout << "lambda_k+1: " << lambda_xp1 << std::endl;*/
		//std::cout  << std::endl;
		output << k << "\t" << error << std::endl;

		lambda_x = lambda_xp1;
		k++;
	}

	output.close();
	return lambda_x;
}

#pragma endregion
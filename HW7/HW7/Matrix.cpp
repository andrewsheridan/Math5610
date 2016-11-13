//Andrew Sheridan
//Math 5610 
//Written in C++
//Matrix.cpp

#include "Matrix.h"
#include "Vector.h"
#include <random>
#include <iostream>
#include <iomanip>


#pragma region Constructors and Initialization
Matrix::Matrix(unsigned size) : rows(size), columns(size)
{
	entries = new Vector[size];
	for (unsigned i = 0; i < size; i++) {
		entries[i] = Vector(size);
		for (unsigned j = 0; j < size; j++) {
			entries[i][j] = 0;
		}
	}
}


Matrix::Matrix(unsigned rowCount, unsigned columnCount) : rows(rowCount), columns(columnCount) {
	//entries = new double*[rowCount];
	entries = new Vector[rowCount];
	for (unsigned i = 0; i < rowCount; i++) {
		//entries[i] = new double[columnCount];
		entries[i] = Vector(columnCount);
		for (unsigned j = 0; j < columnCount; j++) {
			entries[i][j] = 0;
		}
	}
}

///Copy Constructor
Matrix::Matrix(const Matrix &m) : rows(m.rows), columns(m.columns) {
	//entries = new double*[rowsgram];
	entries = new Vector[rows];
	for (unsigned int i = 0; i < rows; i++) {
		//entries[i] = new double[columns];
		entries[i] = Vector(columns);
		for (unsigned int j = 0; j < columns; j++) {
			entries[i][j] = m.entries[i][j];
		}
		//std::cout << i << ": ";
		//entries[i].Print();
	}
}

//Deconstructor
Matrix::~Matrix() {

}

Matrix Matrix::operator=(const Matrix& m) {
	this->rows = m.rows;
	this->columns = m.columns;
	this->entries = new Vector[m.rows];
	for (int i = 0; i < m.rows; i++) {
		this->entries[i] = m.entries[i];
	}

	return *this;
}


///Initializes the matrix to have ones on the main diagonal
void Matrix::InitializeIdentityMatrix() {
	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < i; j++) {
			entries[i][j] = 0;
		}
		entries[i][i] = 1;
		for (unsigned j = i + 1; j < rows; j++) {
			entries[i][j] = 0;
		}
	}
}

///Sets the values of all entries between 0 and 1
void Matrix::InitializeRandom() {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	for (unsigned i = 0; i < columns; i++) {
		for (unsigned j = 0; j < rows; j++) {
			entries[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}
}

///Sets the values of all entries between minValue and maxValue
void Matrix::InitializeRange(double minValue, double maxValue) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(minValue, maxValue); //Desired distribution

	for (unsigned i = 0; i < columns; i++) {
		for (unsigned j = 0; j < rows; j++) {
			entries[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}
}

/// Sets the values of the entries to be diagonally dominant
void Matrix::InitializeDiagonallyDominant() {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < columns; j++) {
			entries[i][j] = dis(generator); //Assign each entry in matrix to random number between 0 and 1
		}
	}

	for (unsigned k = 0; k < rows; k++) {
		entries[k][k] += 10 * rows; //Add 10*n to all diagonal entries
	}
}

#pragma endregion

#pragma region Printing
///Outputs matrix to the console
void Matrix::Print() {
	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < columns; j++) {
			std::cout << std::setw(13) << std::left << entries[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

///Outputs an augmented coefficient matrix to the console
void Matrix::PrintAugmented(Vector v) {
	if (v.GetSize()	 != rows) {
		std::cout << "Vector and Matrix do not match sizes" << std::endl;
		return;
	}

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < columns; j++) {
			std::cout << std::setw(13) << std::left << entries[i][j];
		}
		std::cout << "| " << v[i] << std::endl;
	}
	std::cout << std::endl;
}
#pragma endregion

#pragma region Comparison Operations
bool operator == (const Matrix& A, const Matrix& B) {
	if (A.columns != B.columns) return false;
	if (A.rows != B.rows) return false;

	for (unsigned i = 0; i < A.rows; i++) {
		for (unsigned j = 0; j < A.columns; j++) {
			if (A.entries[i][j] != B.entries[i][j]) {
				return false;
			}
		}
	}
	return true;
}

bool operator != (const Matrix& A, const Matrix& B) {
	if (A.columns != B.columns) return true;
	if (A.rows != B.rows) return true;

	for (unsigned i = 0; i < A.rows; i++) {
		for (unsigned j = 0; j < A.columns; j++) {
			if (A.entries[i][j] != B.entries[i][j]) {
				return true;
			}
		}
	}
	return false;
}

///Checks to see if matrix is symmetric
bool Matrix::IsSymmetric() {
	for (unsigned int i = 0; i < rows; i++) {
		for (unsigned int j = 0; j <= i; j++) {
			if (entries[i][j] != entries[j][i])
				return false;
		}
	}
	return true;
}

#pragma endregion


///Multiplies an nxn matrix A by the vector x
Vector operator *(const Matrix& A, Vector& x) {
	if (A.columns != x.GetSize()) return NULL;

	Vector result(A.rows);
	for (unsigned i = 0; i < A.rows; i++) {
		result[i] = 0;
		for (unsigned j = 0; j < A.columns; j++) {
			result[i] += A.entries[i][j] * x[j];
		}
	}
	return result;
}

///Divides an nxn matrix A by the vector x
Vector operator /(const Matrix& A, Vector& x) {
	if (A.columns != x.GetSize()) return NULL;

	Vector result(A.rows);
	for (unsigned i = 0; i < A.rows; i++) {
		result[i] = 0;
		for (unsigned j = 0; j < A.columns; j++) {
			result[i] += A.entries[i][j] / x[j];
		}
	}
	return result;
}

#pragma endregion

#pragma region operations
/// Returns the transpose of the n by n matrix A
Matrix Matrix::Transpose() {
	Matrix matrix(columns, rows);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < columns; j++) {
			matrix[j][i] = entries[i][j];
		}
	}
	return matrix;
}

///Multiplies an matrix A by matrix B
Matrix operator* (Matrix A, Matrix B) {
	if (A.columns != B.rows) throw "Incompatible sizes";
	
	Matrix matrix(A.rows, B.columns);
	Matrix bTranspose = B.Transpose();
	
	//std::cout << "Starting multiplication. A:" << std::endl;
	/*A.Print();
	std::cout << "B: " << std::endl;
	B.Print();
	std::cout << "Bt: " << std::endl;
	bTranspose.Print();*/

	for (unsigned i = 0; i < A.rows; i++) {
		for (unsigned j = 0; j < B.columns; j++) {
			//matrix[i][j] = DotProduct(A[i], bTranspose[j], A.columns);
			matrix[i][j] = A[i] * bTranspose[j];
		}
	}
	return matrix;
}

///Multiplies an matrix A by matrix B
Matrix operator- (Matrix A, Matrix B) {
	if (A.columns != B.columns && A.rows != B.rows) throw "Incompatible sizes";
	Matrix matrix(A.rows, B.columns);
	for (unsigned i = 0; i < A.rows; i++) {
		for (unsigned j = 0; j < A.columns; j++) {
			matrix[i][j] = A[i][j] - B[i][j];
		}
	}
	return matrix;
}

///Computes the 1-norm of the matrix
double Matrix::OneNorm() {
	double columnMax = 0;
	for (unsigned i = 0; i < rows; i++) {
		double columnSum = 0;
		for (unsigned j = 0; j < columns; j++) {
			columnSum += std::abs(entries[i][j]);
		}
		if (columnSum > columnMax)
			columnMax = columnSum;
	}
	return columnMax;
}

///Computes the inifinity norm of an n by n matrix A
double Matrix::InfinityNorm() {
	double rowMax = 0;
	for (unsigned j = 0; j < columns; j++) {
		double rowSum = 0;
		for (unsigned i = 0; i < rows; i++) {
			rowSum += std::abs(entries[i][j]);
		}
		if (rowSum > rowMax)
			rowMax = rowSum;
	}
	return rowMax;
}



#pragma endregion


#pragma endregion
//Andrew Sheridan
//Math 5610 
//Written in C++
//Matrix.h

#pragma once
#include "Vector.h"

class Matrix{
public: 
	//Initialization and Deconstruction
	Matrix() = default;
	Matrix(unsigned size);
	Matrix(unsigned rowCount, unsigned columnCount);
	Matrix(const Matrix &m);
	Matrix operator = (const Matrix& m);
	~Matrix();
	
	void InitializeIdentityMatrix();
	void InitializeRandom();
	void InitializeRange(double minValue, double maxValue);
	void InitializeDiagonallyDominant();
	
	//Overloaded Operators
	Vector &operator[] (unsigned row) { return entries[row]; }
	friend bool operator == (const Matrix& A, const Matrix& B);
	friend bool operator != (const Matrix& A, const Matrix& B);
	friend Vector operator * (const Matrix& A, Vector& x);
	friend Matrix operator * (Matrix A, Matrix B);
	friend Matrix operator - (Matrix A, Matrix B);

	//Basic Algorithms
	bool IsSymmetric();
	Matrix Transpose();
	double OneNorm();
	double InfinityNorm();

	//Getters and Setters
	unsigned GetRows() { return rows; }
	unsigned GetColumns() { return columns; }
	void SetRows(unsigned r) { rows = r; }
	void SetColumns(unsigned c) { columns = c; }

	//Output
	void Print();
	void PrintAugmented(Vector v);

private:
	Vector* entries; //The entries of the matrix

	unsigned rows;
	unsigned columns;
};

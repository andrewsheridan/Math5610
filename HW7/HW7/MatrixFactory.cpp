//Andrew Sheridan
//Math 5610 
//Written in C++
//MatrixFactory.cpp

#include "MatrixFactory.h"
#include <random>

MatrixFactory* MatrixFactory::m_instance = nullptr;

///Returns the instance of the MatrixFactory singleton
MatrixFactory* MatrixFactory::Instance() {
	if (!m_instance)
		m_instance = new MatrixFactory();

	return m_instance;
}

///Creates an identity matrix
Matrix MatrixFactory::Identity(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	m.InitializeIdentityMatrix();
	return m;
}

///Creates a matrix where every entry is a 1
Matrix MatrixFactory::Ones(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			m[i][j] = 1;
		}
	}
	return m;
}

///Creates a matrix where every entry is a value between 0 and 1
Matrix MatrixFactory::Random(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	m.InitializeRandom();
	return m;
}

///Creates a matrix where every entry is a value between minValue and maxValue
Matrix MatrixFactory::RandomRange(unsigned rows, unsigned columns, double minValue, double maxValue) {
	Matrix m(rows, columns);
	m.InitializeRange(minValue, maxValue);
	return m;
}

///Creates a Diagonally Dominant matrix
Matrix MatrixFactory::DiagonallyDominant(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	m.InitializeDiagonallyDominant();
	return m;
}

///Creates a symmetric matrix, where all entries have values between 0 and 1
Matrix MatrixFactory::Symmetric(unsigned size) {
	Matrix m(size);
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = i; j < size; j++) {
			double value = dis(generator); //Assign each entry in matrix to random number between 0 and 1
			m[i][j] = value;
			m[j][i] = value;
		}
	}

	return m;
}

///Creates a symmetric positive definite matrix
Matrix MatrixFactory::SPD(unsigned size) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	Matrix matrix(size);

	for (unsigned i = 0; i < size; i++) {
		for (unsigned j = i; j < size; j++) {
			double value = dis(generator); //Assign each entry in matrix to random number between 0 and 1
			matrix[i][j] = value;
			matrix[j][i] = value;
		}
	}

	for (unsigned k = 0; k < size; k++) {
		matrix[k][k] += 10 * size; //Add 10*n to all diagonal entries
	}

	return matrix;
}
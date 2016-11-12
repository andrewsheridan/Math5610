//Andrew Sheridan
//Math 5610 
//Written in C++
//MatrixFactory.cpp

#include "MatrixFactory.h"

Matrix MatrixFactory::Identity(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	m.InitializeIdentityMatrix();
	return m;
}

Matrix MatrixFactory::Random(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	m.InitializeRandom();
	return m;
}

Matrix MatrixFactory::RandomRange(unsigned rows, unsigned columns, double minValue, double maxValue) {
	Matrix m(rows, columns);
	m.InitializeRange(minValue, maxValue);
	return m;
}

Matrix MatrixFactory::DiagonallyDominant(unsigned rows, unsigned columns) {
	Matrix m(rows, columns);
	m.InitializeDiagonallyDominant();
	return m;
}

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
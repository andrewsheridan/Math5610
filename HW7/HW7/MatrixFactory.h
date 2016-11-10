#pragma once
#include "Matrix.h"

class MatrixFactory {
public:
	Matrix Identity(unsigned rows, unsigned columns);
	Matrix Random(unsigned rows, unsigned columns);
	Matrix RandomRange(unsigned rows, unsigned columns, double min, double max);
	Matrix DiagonallyDominant(unsigned rows, unsigned columns);
	Matrix Symmetric(unsigned size);
	Matrix SPD(unsigned size);
};
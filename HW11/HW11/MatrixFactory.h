//Andrew Sheridan
//Math 5610 
//Written in C++
//MatrixFactory.h
#pragma once
#include "Matrix.h"

///A singleton which creates new matrices
class MatrixFactory {
public:
	static MatrixFactory* Instance();
	Matrix Identity(unsigned rows, unsigned columns);
	Matrix Ones(unsigned rows, unsigned columns);
	Matrix UpperTriangular(unsigned rows, unsigned columns);
	Matrix LowerTriangular(unsigned rows, unsigned columns);
	Matrix Random(unsigned rows, unsigned columns);
	Matrix RandomRange(unsigned rows, unsigned columns, double min, double max);
	Matrix DiagonallyDominant(unsigned rows, unsigned columns);
	Matrix Symmetric(unsigned size);
	Matrix SPD(unsigned size);
	
private: 
	MatrixFactory() {};
	MatrixFactory(MatrixFactory const&) {};
	MatrixFactory& operator=(MatrixFactory const&) {};
	static MatrixFactory* m_instance;
};
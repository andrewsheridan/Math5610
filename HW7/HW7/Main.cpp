//Andrew Sheridan
//Math 5610 
//Written in C++
//main.cpp

#include "Matrix.h"
#include "Vector.h"
#include "MatrixOperations.h"
#include "MatrixFactory.h"
#include <iostream>

int main() {
	const unsigned size = 3;
	MatrixFactory MF;
	
	//Problem 1
	Matrix matrix1(5, 3);
	
	matrix1[0][0] = 1;
	matrix1[0][1] = 0;
	matrix1[0][2] = 1;
	matrix1[1][0] = 2;
	matrix1[1][1] = 3;
	matrix1[1][2] = 5;
	matrix1[2][0] = 5;
	matrix1[2][1] = 3;
	matrix1[2][2] = -2;
	matrix1[3][0] = 3;
	matrix1[3][1] = 5;
	matrix1[3][2] = 4;
	matrix1[4][0] = -1;
	matrix1[4][1] = 6;
	matrix1[4][2] = 3;
	Vector blah = Vector(5);	
	blah[0] = 4;
	blah[1] = -2;
	blah[2] = 5;
	blah[3] = -2;
	blah[4] = 1;
	matrix1.PrintAugmented(blah);

	Vector result1 = LeastSquares(matrix1, blah);

	//Problem 2 and 3
	//Matrix matrix2(2);
	//matrix2[0][0] = 4;
	//matrix2[0][1] = -2;
	//matrix2[1][0] = 3;
	//matrix2[1][1] = 1;
	//Matrix* QR = GramSchmidt(matrix2);
	//Matrix Q = QR[0];
	//Matrix R = QR[1];
	//matrix2.Print();
	//Q.Print();
	//R.Print();

	//Matrix Product = Q * R;
	//Product.Print();

	//Matrix matrix2b(3);
	//matrix2b[0][0] = 1;
	//matrix2b[0][1] = 1;
	//matrix2b[0][2] = 0;
	//matrix2b[1][0] = 1;
	//matrix2b[1][1] = 0;
	//matrix2b[1][2] = 1;
	//matrix2b[2][0] = 0;
	//matrix2b[2][1] = 1;
	//matrix2b[2][2] = 1;

	//Matrix* QRb = GramSchmidt(matrix2b);
	//Matrix Qb = QRb[0];
	//Matrix Rb = QRb[1];
	//Qb.Print();
	//Rb.Print();




	int input;
	std::cin >> input;

	return 0;
}
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
	//Problem 1
	/*Matrix matrix1(5, 3);
	
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
	Vector vector1 = Vector(5);	
	vector1[0] = 4;
	vector1[1] = -2;
	vector1[2] = 5;
	vector1[3] = -2;
	vector1[4] = 1;

	Vector result1 = LeastSquares(matrix1, vector1);
	result1.Print();*/

	//Problem 2 
	Matrix matrix2 = MatrixFactory::Instance()->DiagonallyDominant(5, 5);
	Matrix* QR = GramSchmidt(matrix2);
	Matrix Q = QR[0];
	Matrix R = QR[1];
	
	std::cout << "Starting Matrix: " << std::endl;
	matrix2.Print();
	std::cout << "Q: " << std::endl;
	Q.Print();
	std::cout << "R: " << std::endl;
	R.Print();
	
	std::cout << "Q * R" << std::endl;
	Matrix TestResult = Q * R;
	TestResult.Print();

	Matrix Norm = matrix2 - TestResult;
	double normValue = Norm.OneNorm();
	std::cout << "||A - (Q * R)|| : " << normValue << std::endl;
	//Problem 3
	/*Matrix matrix3(2);
	matrix3[0][0] = 4;
	matrix3[0][1] = -2;
	matrix3[1][0] = 3;
	matrix3[1][1] = 1;
	Matrix* QR = GramSchmidt(matrix3);
	Matrix Q = QR[0];
	Matrix R = QR[1];
	matrix3.Print();
	Q.Print();
	R.Print();

	Matrix Product = Q * R;
	Product.Print();

	Matrix matrix3b(3);
	matrix3b[0][0] = 1;
	matrix3b[0][1] = 1;
	matrix3b[0][2] = 0;
	matrix3b[1][0] = 1;
	matrix3b[1][1] = 0;
	matrix3b[1][2] = 1;
	matrix3b[2][0] = 0;
	matrix3b[2][1] = 1;
	matrix3b[2][2] = 1;

	Matrix* QRb = GramSchmidt(matrix3b);
	Matrix Qb = QRb[0];
	Matrix Rb = QRb[1];
	Qb.Print();
	Rb.Print();*/


	int input;
	std::cin >> input;

	return 0;
}
//Andrew Sheridan
//Math 5610 
//Written in C++
//main.cpp

#include "Matrix.h"
#include "Vector.h"
#include "MatrixOperations.h"
#include "MatrixFactory.h"
#include <iostream>
#include <fstream>
#include <iomanip>

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
	/*Matrix matrix2 = MatrixFactory::Instance()->DiagonallyDominant(5, 5);
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
	std::cout << "||A - (Q * R)|| : " << normValue << std::endl;*/

	//Problem 3
	//Matrix matrix3(2);
	/*matrix3[0][0] = 4;
	matrix3[0][1] = -2;
	matrix3[1][0] = 3;
	matrix3[1][1] = 1;
	Matrix* QR = GramSchmidt(matrix3);
	Matrix Q = QR[0];
	Matrix R = QR[1];
	std::cout << "First Matrix" << std::endl;
	matrix3.Print();
	std::cout << "Q:" << std::endl;
	Q.Print();
	std::cout << "R:" << std::endl;
	R.Print();

	Matrix Product = Q * R;
	std::cout << "QR" << std::endl;
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

	std::cout << "Matrix two" << std::endl;
	matrix3b.Print();
	std::cout << "Q:" << std::endl;
	Qb.Print();
	std::cout << "R:" << std::endl;
	Rb.Print();
	Matrix secondProduct = Qb*Rb;
	std::cout << "QR:" << std::endl;
	secondProduct.Print();*/

	///Problem 4
	//const int size4 = 5;
	////Creates a diagonally dominant matrix
	//Matrix matrix4 = MatrixFactory::Instance()->DiagonallyDominant(size4, size4);
	//Vector vector4(size4);
	//vector4.InitializeRandomEntries();
	//Vector vector4a = vector4;
	//Vector vector4b = vector4;
	//Matrix matrix4a = matrix4;
	//Matrix matrix4b = matrix4;

	////The first entry in this array is Q, the second, R
	//Matrix* QR4 = GramSchmidt(matrix4a);
	//Vector QTy = QR4[0].Transpose() * vector4a;
	//Vector x = BackSubstitution(QR4[1], QTy);

	//GaussianElimination(matrix4b, vector4b);
	//Vector gaussianEliminationResults = BackSubstitution(matrix4b, vector4b);

	//std::cout << "Test System: " << std::endl;
	//matrix4.PrintAugmented(vector4);

	//std::cout << "Results of gaussian elimination and back substitution: " << std::endl;
	//gaussianEliminationResults.Print();

	//std::cout << "Results of QR factorization and back substitution: " << std::endl;
	//x.Print();

	//Vector difference = gaussianEliminationResults - x;
	//std::cout << "L2 norm of the difference:" << difference.L2Norm() << std::endl;

	///Problem 5
	//for (int size5 = 10; size5 <= 160; size5 *= 2) {
	//	//Creates a diagonally dominant matrix
	//	Matrix matrix5 = MatrixFactory::Instance()->DiagonallyDominant(size5, size5);
	//	Vector onesVector(size5);
	//	Vector vector5(size5);
	//	onesVector.InitializeAllOnes();
	//	vector5 = matrix5 * onesVector;

	//	//The first entry in this array is Q, the second, R
	//	Matrix* QR5 = GramSchmidt(matrix5);
	//	Vector QTy = QR5[0].Transpose() * vector5;
	//	Vector Result = BackSubstitution(QR5[1], QTy);

	//	/*std::cout << "Results with size " << size5 << std::endl;
	//	Result.Print();*/

	//	Vector difference = Result - onesVector;
	//	std::cout << "Current size: " << size5 << std::endl;
	//	std::cout << "L2 Error : " << difference.L2Norm() << std::endl;
	//}
	
	///Problem 6
	Vector v6(6);
	Vector b6(6);
	for (int i = 0; i < 6; i++) {
		v6[i] = i * 0.2;
		b6[i] = std::rand();
	}
	LSFit(v6, b6, 4);

	///Problem 7
	std::ofstream output("output.txt");

	Vector t(11);
	Vector b(11);
	t[0] = 0.0;	
	t[1] = 0.1;
	t[2] = 0.2;
	t[3] = 0.3;
	t[4] = 0.4;
	t[5] = 0.5;
	t[6] = 0.6;
	t[7] = 0.7;
	t[8] = 0.8;
	t[9] = 0.9;
	t[10] = 1;

	b[0] = 0.9;
	b[1] = 1.01;
	b[2] = 1.05;
	b[3] = .97;
	b[4] = .98;
	b[5] = .95;
	b[6] = 0.01;
	b[7] = -.1;
	b[8] = .02;
	b[9] = -.1;
	b[10] = 0;

	const int N = 5;
	Vector* coefs = new Vector[N];
	for (int i = 0; i <= N; i++) {
		coefs[i] = LSFit(t, b, i);
	}
	int size = 1 / .01;
	Vector T(size);
	for (int i = 0; i < size; i++) {
		T[i] = .01 * i;
	}

	Matrix Z = MatrixFactory::Instance()->Ones(N, size);
	
	for (int n = 0; n < N; n++) {
		for (int i = 0; i < size; i++) {
			Z[n] = Z[n] * coefs[n][n];
			for (int j = n - 2; j >= 0; j--) {
				Z[n][j] = (Z[n] * T) + coefs[n][j];
			}
		}
	}

	output << std::setw(13) << "t" << std::setw(13) << "p1" << std::setw(13) << "p2" 
		<< std::setw(13) << "p1" << std::setw(13) << "p3" << std::setw(13) << "p4" 
		<< std::setw(13) << "p5" << std::endl;
	for (int j = 0; j < size; j++) {
		output << std::setw(13) << T[j];
		for (int n = 0; n < N; n++) {
			output << std::setw(13) << Z[n][j];
		}
		output << std::endl;
	}
	output.close();

	int input;
	std::cin >> input;

	return 0;
}
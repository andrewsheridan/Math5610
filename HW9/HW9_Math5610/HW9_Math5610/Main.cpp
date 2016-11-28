//Andrew Sheridan
//Math 5610 
//Written in C++
//main.cpp

#include "../../../HW7/HW7/Matrix.h"
#include "../../../HW7/HW7/MatrixFactory.h"
#include "../../../HW7/HW7/MatrixOperations.h"
#include "../../../HW7/HW7/Vector.h"

int main() {
	//Problem 1 : Jacobi Iteration
	//Matrix m1(3);
	//m1[0][0] = 7;
	//m1[0][1] = 3;
	//m1[0][2] = 1;
	//m1[1][0] = -3;
	//m1[1][1] = 10;
	//m1[1][2] = 2;
	//m1[2][0] = 1;
	//m1[2][1] = 7;
	//m1[2][2] = -15;

	//Vector b1(3);
	//b1[0] = 3;
	//b1[1] = 4;
	//b1[2] = 2;

	//Vector x1(3);

	//int maxIter = 10000;
	//double tolerance = 0.00001;
	//std::cout << "The example matrix given on page 175 of the textbook: " << std::endl;
	//m1.PrintAugmented(b1);

	//std::cout << "The results of Jacobi Iteration with max iterations of " << maxIter << " and tolerance of" << tolerance << std::endl;
	//Vector result1 = JacobiIteration(m1, x1, b1, maxIter, tolerance);
	//result1.Print();

	//std::cout << "The results of Gaussian Elimination and Back Substitution" << std::endl;
	//GaussianElimination(m1, b1);
	//Vector actual1 = BackSubstitution(m1, b1);
	//actual1.Print();

	////Problem 2 : Gauss-Seidel
	//Matrix m2(3);
	//m2[0][0] = 7;
	//m2[0][1] = 3;
	//m2[0][2] = 1;
	//m2[1][0] = -3;
	//m2[1][1] = 10;
	//m2[1][2] = 2;
	//m2[2][0] = 1;
	//m2[2][1] = 7;
	//m2[2][2] = -15;

	//Vector b2(3);
	//b2[0] = 3;
	//b2[1] = 4;
	//b2[2] = 2;

	//Vector x2(3);
	//std::cout << "The example matrix given on page 175 of the textbook: " << std::endl;
	//m2.PrintAugmented(b2);

	//std::cout << "The results of Gauss Seidel with max iterations of " << maxIter << " and tolerance of" << tolerance << std::endl;
	//Vector result2 = GaussSeidel(m2, x2, b2, maxIter, tolerance);
	//result2.Print();

	//std::cout << "The results of Gaussian Elimination and Back Substitution" << std::endl;
	//GaussianElimination(m2, b2);
	//Vector actual2 = BackSubstitution(m2, b2);
	//actual2.Print();

	////Problem 3: Graph Differences
	//int maxIter3 = 10000;
	//double tol3 = 0.00001;
	//for (int i = 4; i <= 160; i *= 2) {
	//	Matrix m3 = MatrixFactory::Instance()->DiagonallyDominant(i, i);
	//	Vector onesVector(i);
	//	onesVector.InitializeAllOnes();
	//	Vector v3 = m3 * onesVector;

	//	Vector zeroes(i);
	//	Vector result3 = JacobiIteration(m3, zeroes, v3, maxIter3, tol3);
	//	std::cout << "Result of Jacobi iteration on DD matrix size " << i << std::endl;
	//	result3.Print();

	//	Vector result3b = GaussSeidel(m3, zeroes, v3, maxIter3, tol3);
	//	std::cout << "Result of Gauss Seidel on DD matrix size " << i << std::endl;
	//	result3b.Print();
	//}

	////Problem 4: LU Factorization vs Iterative Methods
	//int maxIter4 = 10000;
	//double tol4 = 0.000000001;
	//for (int i = 5; i <= 160; i *= 2) {
	//	Matrix m4 = MatrixFactory::Instance()->DiagonallyDominant(i, i);
	//	Vector onesVector(i);
	//	onesVector.InitializeAllOnes();
	//	Vector v4 = m4 * onesVector;

	//	Vector zeroes(i);
	//	Vector result4 = JacobiIteration(m4, zeroes, v4, maxIter4, tol4);

	//	Vector result4b = GaussSeidel(m4, zeroes, v4, maxIter4, tol4);

	//	LUFactorization(m4, v4);
	//	Vector result4c = BackSubstitution(m4, v4);
	//}

	/*for (int i = 5; i <= 160; i *= 2) {
		Matrix m4 = MatrixFactory::Instance()->Random(i, i);
		Vector onesVector(i);
		onesVector.InitializeAllOnes();
		Vector v4 = m4 * onesVector;

		Vector zeroes(i);
		Vector result4 = JacobiIteration(m4, zeroes, v4, maxIter4, tol4);
		std::cout << "Result of Jacobi iteration on DD matrix size " << i << std::endl;
		result4.Print();

		Vector result4b = GaussSeidel(m4, zeroes, v4, maxIter4, tol4);
		std::cout << "Result of Gauss Seidel on DD matrix size " << i << std::endl;
		result4b.Print();

		LUFactorization(m4, v4);
		Vector result4c = BackSubstitution(m4, v4);
		std::cout << "Result of LU Factorization and Back Substitution on DD matrix size " << i << std::endl;
		result4c.Print();
	}*/

	//Problem 5 : Conjugate Gradient Method
	Matrix m5(3);
	m5[0][0] = 7;
	m5[0][1] = 3;
	m5[0][2] = 1;
	m5[1][0] = 3;
	m5[1][1] = 10;
	m5[1][2] = 2;
	m5[2][0] = 1;
	m5[2][1] = 2;
	m5[2][2] = 15;

	Vector b5(3);
	b5[0] = 28;
	b5[1] = 31;
	b5[2] = 22;

	Vector x5(3);

	std::cout << "Matrix and RHS given on page 184." << std::endl;
	m5.PrintAugmented(b5);
	Vector result5 = ConjugateGradient(m5, b5, x5, 0.000001);
	std::cout << "Result of conjugate gradient method: " << std::endl;
	result5.Print();

	for (int i = 5; i <= 160; i *= 2) {
		Matrix m5 = MatrixFactory::Instance()->SPD(i);
		Vector ones(i);
		ones.InitializeAllOnes();
		Vector v5 = m5 * ones;
		Vector zeroes(i);

		Vector result5 = ConjugateGradient(m5, v5, zeroes, 0.0000001);
		result5.Print();
	}
	

	int input;
	std::cin >> input;
	return 0;
}
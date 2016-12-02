#include "../../HW7/HW7/Matrix.h";
#include "../../HW7/HW7/MatrixFactory.h";
#include "../../HW7/HW7/MatrixOperations.h";
#include "../../HW7/HW7/Vector.h";
#include "Error.h";

int main() {
	int size = 4;

	Matrix matrix = MatrixFactory::Instance()->Random(size, size);
	Vector vector(size);
	vector.InitializeAllOnes();
	vector = matrix * vector;
	Vector resultVector = vector;
	Matrix reducedMatrix = GaussianEliminationWithScaledPivoting(matrix, resultVector);

	std::cout << "Gaussian Elimination With Scaled Pivoting" << std::endl;
	std::cout << "Test Matrix" << std::endl;
	matrix.Print();
	std::cout << "Test vector " << std::endl;
	vector.Print();
	std::cout << "Result of gaussian elimination " << std::endl;
	reducedMatrix.PrintAugmented(resultVector);

	//problem 1: Power Method
	Matrix A(2);
	A[0][0] = 2;
	A[0][1] = -12;
	A[1][0] = 1;
	A[1][1] = -5;

	Vector x_0(2);
	x_0[0] = 1;
	x_0[1] = 1;
	std::cout << "Test Matrix: " << std::endl;
	A.PrintAugmented(x_0);
	double lambda = PowerMethod(A, x_0, 0.0001, 1000);
	std::cout << "Lambda: " << lambda << std::endl;
	double relativeError = realRelative(-2.0, lambda);
	std::cout << "Relative Error: " << relativeError << std::endl;
	double absoluteError = realAbsolute(-2.0, lambda);
	std::cout << "Absolute Error: " << absoluteError << std::endl;

	for (int i = 5; i <= 160; i *= 2) {
		Matrix m1 = MatrixFactory::Instance()->DiagonallyDominant(i, i);
		Vector onesVector(i);
		onesVector.InitializeAllOnes();
		Vector v1 = m1 * onesVector;

		Vector zeroes(i);

		double result1a = PowerMethod(m1, onesVector, 0.000001, 10000);
		std::cout << "Result of Power Method on DD matrix size " << i << std::endl;
		std::cout << result1a << std::endl;
	}

	//problem 3: Inverse Power Method
	Matrix A3(2);
	A3[0][0] = 2;
	A3[0][1] = -12;
	A3[1][0] = 1;
	A3[1][1] = -5;

	Vector x3_0(2);
	x3_0[0] = 1;
	x3_0[1] = 1;

	double lambda3 = InversePowerMethod(A3, x3_0, 0.000001, 1000);
	std::cout << "Lambda: " << lambda3 << std::endl;
	double relativeError3 = realRelative(-1.0, lambda3);
	std::cout << "Relative Error: " << relativeError3 << std::endl;
	double absoluteError3 = realAbsolute(-1.0, lambda3);
	std::cout << "Absolute Error: " << absoluteError3 << std::endl;

	//Problem 4 : Inverse Power Method Analysis
	for (int i = 5; i <= 80; i *= 2) {
		Matrix m4 = MatrixFactory::Instance()->DiagonallyDominant(i, i);
		Vector onesVector(i);
		onesVector.InitializeAllOnes();
		Vector v4 = m4 * onesVector;

		Vector zeroes(i);
		std::cout << "Test matrix: " << std::endl;
		double result4a = InversePowerMethod(m4, onesVector, 0.000000001, 10000);
		std::cout << "Result of Inverse Power Method on DD matrix size " << i << std::endl;
		std::cout << result4a << std::endl;
	}

	int input = 0;
	std::cin >> input;

	return 0;
}
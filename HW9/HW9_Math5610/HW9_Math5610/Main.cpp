#include "../../../HW7/HW7/Matrix.h"
#include "../../../HW7/HW7/MatrixFactory.h"
#include "../../../HW7/HW7/MatrixOperations.h"
#include "../../../HW7/HW7/Vector.h"

int main() {
	//Problem 1 : Jacobi Iteration
	Matrix m1(3);
	m1[0][0] = 7;
	m1[0][1] = 3;
	m1[0][2] = 1;
	m1[1][0] = -3;
	m1[1][1] = 10;
	m1[1][2] = 2;
	m1[2][0] = 1;
	m1[2][1] = 7;
	m1[2][2] = -15;

	Vector b1(3);
	b1[0] = 3;
	b1[1] = 4;
	b1[2] = 2;

	Vector x1(3);

	std::cout << "The results of Jacobi Iteration with max iterations of 100 and tolerance of 0.0001." << std::endl;
	Vector result1 = JacobiIteration(m1, x1, b1, 100, 0.0001);
	result1.Print();

	std::cout << "The results of Gaussian Elimination and Back Substitution" << std::endl;
	GaussianElimination(m1, b1);
	Vector actual1 = BackSubstitution(m1, b1);
	actual1.Print();







	int input;
	std::cin >> input;
	return 0;
}
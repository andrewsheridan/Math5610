#include "../../HW7/HW7/Matrix.h";
#include "../../HW7/HW7/MatrixFactory.h";
#include "../../HW7/HW7/MatrixOperations.h";
#include "../../HW7/HW7/Vector.h";

int main() {
	//problem 1: Power Method
	Matrix A(2);
	A[0][0] = 2;
	A[0][1] = -12;
	A[1][0] = 1;
	A[1][1] = -5;

	Vector x_0(2);
	x_0[0] = 1;
	x_0[1] = 1;

	double lambda = PowerMethod(A, x_0, 0.0001, 1000);
	std::cout << "Lambda: " << lambda << std::endl;

	int input = 0;
	std::cin >> input;

	return 0;
}
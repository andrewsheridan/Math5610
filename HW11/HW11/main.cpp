// Andrew Sheridan
//Math 5610 
//Written in C++
//main.cpp
#include "Matrix.h"
#include "MatrixFactory.h"
#include "MatrixOperations.h"
#include "Vector.h"
#include "Error.h"
#include "Complex.h"
#include "ApproximationMethods.h"

double ft(double t) {
	return (3 * t*t) - (2 * t*t*t);
}

double ft2(double t) {
	return sin(3.14159265358979 * t );
}

int main() {
	Vector x(5);
	x[0] = 2.5;
	x[1] = 3;
	x[2] = 4.5;
	x[3] = 4.75;
	x[4] = 6;

	Vector fx(5);
	fx[0] = 8.85;
	fx[1] = 11.45;
	fx[2] = 20.66;
	fx[3] = 22.85;
	fx[4] = 38.60;

	Matrix dd = DividedDifference(x, fx);
	dd.Print();

	for (int n = 4; n <= 32; n *= 2) {
		Vector* XFx = EvaluateFunctionAtPointsInRange(ft, -1, 2, n);
		x = XFx[0];
		fx = XFx[1];
		x.Print();
		fx.Print();

		dd = DividedDifference(x, fx);
		dd.Print();
		std::cout << "Coefficients: ";
		Vector c = GetNewtonFormCoefficients(dd);
		c.Print();
		Vector px(fx.GetSize());
		for (int i = 0; i < x.GetSize(); i++) {
			px[i] = InterpolatingPolynomial(x[i], x, c);
		}
		std::cout << "x: ";
		x.Print();
		std::cout << "f(x): ";
		fx.Print();
		std::cout << "p(x): ";
		px.Print();

	}
	for (int n = 4; n <= 32; n *= 2) {
		Vector* XFx = EvaluateFunctionAtPointsInRange(ft2, 0, 3, n);
		x = XFx[0];
		fx = XFx[1];
		x.Print();
		fx.Print();

		dd = DividedDifference(x, fx);
		dd.Print();
		std::cout << "Coefficients: ";
		Vector c = GetNewtonFormCoefficients(dd);
		c.Print();
		Vector px(fx.GetSize());
		for (int i = 0; i < x.GetSize(); i++) {
			px[i] = InterpolatingPolynomial(x[i], x, c);
		}
		std::cout << "x: ";
		x.Print();
		std::cout << "f(x): ";
		fx.Print();
		std::cout << "p(x): ";
		px.Print();
	}

	int input = 0;
	std::cin >> input;

	return 0;
}

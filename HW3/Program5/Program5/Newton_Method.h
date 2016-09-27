// Andrew Sheridan
// Math 5610
// 09/26/15

//Newton_Method.h
#pragma once
#include <cmath>

double newton(double (*f)(double), double (*df)(double), double x0, double tol, unsigned max_iterations) {
	double xk = 99999999; //x sub k
	double xkp1 = x0; //x sub k+1
	//Loops until we reach the maximum iterations or we are within the input tolerance.
	unsigned k = 0;
	while (k < max_iterations && std::abs(xkp1 - xk) > tol) {
		xk = xkp1;
		xkp1 = xk - (f(xk) / df(xk));
		k++;
	}
	std::cout << "Iterations: " << k << std::endl;
	return xkp1;
}
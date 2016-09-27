#pragma once
#include <cmath>

double newton(double (*f)(double), double (*df)(double), double x0, double tol, unsigned max_iterations) {
	double xk = 99999999; //x sub k
	double xkp1 = x0; //x sub k+1
	//Loops until we reach the maximum iterations or we are within the input tolerance.
	for (unsigned k = 0; k < max_iterations && std::abs(xkp1 - xk) > tol; k++) {
		xk = xkp1;
		xkp1 = xk - (f(xk) / df(xk));
	}

	return xkp1;
}
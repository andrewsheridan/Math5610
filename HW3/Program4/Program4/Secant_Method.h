#pragma once
#include<cmath>

// Secant_Method.h
double secant_method(double (*f)(double), double x0, double x1, double tol, unsigned max_iterations) {
	double xkm1; //x sub k-1
	double xk = x0; //x sub k
	double xkp1 = x1; //x sub k+1
	
	//Loops until we reach the maximum iterations or we are within the input tolerance.
	for (unsigned k = 0; k < max_iterations && std::abs(xkp1 - xk) > tol; k++) {
		xkm1 = xk;
		xk = xkp1;
		xkp1 = xk - (f(xk) * (xk - xkm1)) / (f(xk) - f(xkm1));
	}

	return xkp1;
}

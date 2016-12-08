#pragma once
// Andrew Sheridan
// Math 5610
// 12/8/16

#include <iostream>
#include <cmath>


//From 3.1
///The Bisection Method
///f: The function f, which takes a double and returns a double
///		Input: double
///		Output: double
///a: Left starting point
///b: Right starting point
///tol: The algorithm's tolerance
///max_iterations: The maximimum number of loops before exit
double bisection(double(*f)(double), double a, double b, double tol, unsigned max_iterations) {
	if (a > b) { // If a is greater than b, switch them.
		double temp = b;
		b = a;
		a = temp;
	}
	double fa = f(a);
	double fb = f(b);
	if (a >= b || (fa * fb >= 0) || tol <= 0 || max_iterations < 1) // Validate input
	{
		std::cout << "This input is not valid." << std::endl;
		return 0;
	}

	for (int k = 0; k < max_iterations; k++) { //Iterate n times, each iteration bisecting once.
		double p = (a + b) / 2;
		double fp = f(p);
		if (fa * fp < 0) {
			b = p;
			fb = fp;
		}
		else {
			a = p;
			fa = fp;
		}
	}
	return (a + b) / 2;
}

//From 3.2
///The Secant Method
///g: The function g, which takes a double and returns a double
///		Input: double
///		Output: double
///x0: Starting Point
///tol: The algorithm's tolerance
///max_iterations: The maximimum number of loops before exit
double fixed_point(double(*g)(double), double x0, unsigned max_iterations, double tol) {
	// Check input. Other inputs will work properly because of their types. 
	if (tol < 0) {
		std::cout << "The tolerance must be positive." << std::endl;
	}

	double x_i = 99999999;//x sub i
	double x_ip1 = x0; //x sub i+1

					   // computes x of i+1 until the max interations are reached or we are within the tolerance.
	for (int i = 0; i < max_iterations || std::abs(x_ip1 - g(x_i)) > tol; i++) {
		x_i = x_ip1;
		x_ip1 = g(x_i);
	}
	return x_ip1;
}

//From 3.4
///The Secant Method
///f: The function f, which takes a double and returns a double
///		Input: double
///		Output: double
///x0: First Point
///x1: Second Point
///tol: The algorithm's tolerance
///max_iterations: The maximimum number of loops before exit
double secant_method(double(*f)(double), double x0, double x1, double tol, unsigned max_iterations) {
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

//From 3.3
///The Newton Method
///f: The function f, which takes a double and returns a double
///		Input: double
///		Output: double
///df: The derivative of function f
///		Input: double
///		Output: double
///x0: The initial guess
///tol: The algorithm's tolerance
///max_iterations: The maximimum number of loops before exit
double newton_method(double(*f)(double), double(*df)(double), double x0, double tol, unsigned max_iterations) {
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

//From 3.5
///The Hybrid Method, which executes the bisection method until newton method starts to converge
///f: The function f, which takes a double and returns a double
///		Input: double
///		Output: double
///df: The derivative of function f
///		Input: double
///		Output: double
///a: Left bound of initial guess
///b: Right bound of initial guess
///bisection_iterations: The number of bisections executed before attempting newton method
///tol: The algorithm's tolerance
///max_iterations: The maximimum number of loops before exit
double hybrid_method(double(*f)(double), double(*df)(double), double a, double b, unsigned int bisection_iterations, unsigned int max_iterations, double tol) {

	double fa = f(a);
	double fb = f(b);
	// Validate input
	if (a >= b || (fa * fb >= 0) || tol <= 0 || bisection_iterations < 1 || max_iterations < 1)
	{
		std::cout << "This input is not valid." << std::endl;
		return 0;
	}

	unsigned int totalIterationCount = 0;
	unsigned int bisectionIterationCount = 0;
	bool useNewton = false;
	double p;
	while (totalIterationCount < max_iterations && !useNewton && (b - a) > tol) {
		p = (a + b) / 2;
		double fp = f(p);
		if (fa * fp < 0) {
			b = p;
			fb = fp;
		}
		else {
			a = p;
			fa = fp;
		}

		bisectionIterationCount++;
		totalIterationCount++;

		// If we've done the set number of bisections, try newton's method
		if (bisectionIterationCount == bisection_iterations) {
			double newtonResult = newton(f, df, p, tol, 1);
			// If newton's method was more efficient than bisection, start using newton's method. 
			if (std::abs(newtonResult - p) < std::abs(b - a)) {
				useNewton = true;
				p = newtonResult;
				totalIterationCount++;
			}
			//Else, start bisection again.
			else {
				bisectionIterationCount = 0;
			}
		}
	}
	// If we have begun using newton's method, iterate through newton's method until we've reached the max number of iterations or tolerance is met. 
	if (useNewton) {
		double previous = 999999;
		while (totalIterationCount < max_iterations && std::abs(p - previous) > tol) {
			previous = p;
			p = newton(f, df, p, tol, max_iterations - totalIterationCount);
			totalIterationCount++;
		}
	}

	return p;
}
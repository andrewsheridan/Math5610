//Andrew Sheridan
//Math 5610
//09/26/15

#pragma once

#include <math.h>
#include <iostream>

using namespace std;

double bisection(double (*f)(double), double a, double b, double atol) {
	if (a > b) { // If a is greater than b, switch them.
		double temp = b;
		b = a;
		a = temp;
	}
	double fa = f(a);
	double fb = f(b);
	if (a >= b || (fa * fb >= 0) || atol <= 0) // Validate input
	{
		std::cout << "This input is not valid." << std::endl;
		return 0;
	}

	int n = ceil(log2(b - a) - log2(2 * atol)); //Compute number of iterations.
	unsigned k = 0;
	while (k < n) {
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
		k++;
	}
	std::cout << "Iterations: " << k;
	return (a + b) / 2;
}

double bisection_max_iterations(double(*f)(double), double a, double b, double atol, unsigned max_iterations) {
	if (a > b) { // If a is greater than b, switch them.
		double temp = b;
		b = a;
		a = temp;
	}
	double fa = f(a);
	double fb = f(b);
	if (a >= b || (fa * fb >= 0) || atol <= 0 || max_iterations < 1) // Validate input
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

double* bisect_once(double(*f)(double), double a, double b) {
	double fa = f(a);
	double fb = f(b);
	if (a >= b || (fa * fb >= 0)) // Validate input
	{
		std::cout << "This input is not valid." << std::endl;
		return 0;
	}

	double p = (a + b) / 2;
	double fp = f(p);
	if (fa * fp < 0) {
		return new double[a, p];
	}
	else {
		return new double[p, b];
	}
}

double bisect_recursive(double(*f)(double), double a, double b, double fa, double fb, double atol) {
	double p = (a + b) / 2;
	if ((b - a) < atol) {
		return p;
	}
	else {
		double fp = f(p);
		if (fa * fp < 0) {
			b = p;
			fb = fp;
		}
		else {
			a = p; 
			fa = fp;
		}
		p = bisect_recursive(f, a, b, fa, fb, atol);
	}
}
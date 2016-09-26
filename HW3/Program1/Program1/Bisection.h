#pragma once

#include <math.h>

using namespace std;

double bisection(double (*f)(double), double a, double b, double atol) {
	if (a > b) {
		double temp = b;
		b = a;
		a = temp;
	}
	double fa = f(a);
	double fb = f(b);
	if (a >= b || (fa * fb >= 0) || atol <= 0) // Validate input
		throw;

	int n = ceil(log2(b - a) - log2(2 * atol));

	for (int k = 0; k < n; k++) { //Iterate n times, each iteration bisecting once.
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
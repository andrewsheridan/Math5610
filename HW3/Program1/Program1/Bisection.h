#pragma once

#include <cmath>

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

	for (int k = 0; k < n; k++) {
		double p = (a + b) / 2;
		double fp = f(p);
		if (fa * fb < 0) {
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
#pragma once

#include <cmath>

using namespace std;

double bisection(double (*f)(), double a, double b, double fa, double fb, double atol) {
	if (a >= b || (fa * fb <= 0) || atol <= 0)
		throw "Invalid input";

	int n = ceil(log2(b - a) - log2(2 * atol));

	for (int k = 0; k < n; k++) {
		double p = (a + b) / 2;
		fp = f(p);
		
	}
}
#pragma once

double newton(double (*f)(double), double (*df)(double), double x0, double tol, unsigned n) {
	double xk = x0;
	for (unsigned k = 0; k < n; k++) {
		double xkp1 = xk - (f(xk) / df(xk));

	}

}
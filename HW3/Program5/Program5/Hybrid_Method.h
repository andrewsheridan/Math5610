#pragma once
//Hybrid_Method.h
#include <math.h>
#include <iostream>
#include "Newton_Method.h"

double hybrid(double(*f)(double), double(*df)(double), double a, double b, unsigned int bisection_iterations, unsigned int max_iterations, double tol) {
	
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
	while (totalIterationCount < max_iterations && !useNewton && (b-a) > tol) {
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
		if(bisectionIterationCount == bisection_iterations){
			double newtonResult = newton(f, df, p, tol, 1);
			// If newton's method was more efficient than bisection, start using newton's method. 
			if (std::abs(newtonResult - p) < std::abs(b - a)){
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
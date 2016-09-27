// Andrew Sheridan
// Math 5610
// 09/26/15

#pragma once
#include <cmath>
#include <iostream>

//Fixed_Point.h

double fixed_point(double(*g)(double), double x0, unsigned max_iterations, double atol) {

	// Check input. Other inputs will work properly because of their types. 
	if (atol < 0) { 
		std::cout << "The tolerance must be positive." << std::endl;
	}

	double x_i = 99999999;//x sub i
	double x_ip1 = x0; //x sub i+1

	// computes x of i+1 until the max interations are reached or we are within the tolerance.
	for (int i = 0; i < max_iterations || std::abs(x_ip1 - g(x_i)) > atol; i++) {
		x_i = x_ip1;
		x_ip1 = g(x_i);
	}
	return x_ip1;
}
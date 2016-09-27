#pragma once
#include<cmath>

double fixed_point(double(*g)(double), double x0, int max_iterations, double atol) {
	double x_i = 99999999;//x sub i
	double x_ip1 = x0; //x sub i+1
	// computes x of i+1 until the max interations are reached or we are within the tolerance.
	for (int i = 0; i < max_iterations || std::abs(x_ip1 - g(x_i)) > atol; i++) {
		x_i = x_ip1;
		x_ip1 = g(x_i);
	}
	return x_ip1;
}
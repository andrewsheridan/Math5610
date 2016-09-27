//Andrew Sheridan
//Math 5610
//Written in C++

#include <iostream>
#include "Newton_Method.h"
#include <cmath>

//A test function, f(x)
double f(double x) {
	return std::pow(x, 2) - 2;
}

// A test function, f'(x), to be paired with f(x)
double df(double x) {
	return 2 * x;
}

int main(void) {
	double result = newton(f, df, 1, 0.0000001, 10);
	std::cout << "result: " << result << std::endl;
	return 0;
}

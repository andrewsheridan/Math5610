//Andrew Sheridan
//Math 5610
//Written in C++

// Main.cpp
#include <iostream>
#include "Secant_Method.h"
#include <cmath>

//A test function, f(x), as given in Example 3.8
double f(double x) {
	return 2 * std::cosh(x/4) - x;
}

int main(void) {
	double result = secant_method(f, 2, 1, std::pow(10, -8), 10);
	std::cout << "result: " << result << std::endl;
	std::cin >> result;
	return 0;
}

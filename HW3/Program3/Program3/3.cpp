//Andrew Sheridan
//Math 5610
//Written in C++

//Main.cpp
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

// Functions from Example 3.7
double g(double x) {
	return 2 * cosh(x / 4) - x;
}

double dg(double x) {
	return 0.5 * sinh(x / 4) - 1;
}

double pten(double x) {
	return (x - 1)*(x - 1)*std::exp(x);
}

double ptenprime(double x) {
	return std::exp(x) * (x - 1) * (x + 1);
}

int main(void) {
	double result = newton(f, df, 1, 0.0000001, 10);
	std::cout << "result: " << result << std::endl;

	double result2 = newton(g, dg, 2, pow(10, -8), 10);
	std::cout << "result2: " << result2 << std::endl;

	double problem10 = newton(pten, ptenprime, 2, pow(10, -8), 15);
	std::cout << "problem10: " << problem10 << std::endl;
	std::cin >> result;

	return 0;
}

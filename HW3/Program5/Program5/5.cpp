//Main.cpp
#include <iostream>
#include "Hybrid_Method.h"
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
	double result = hybrid(f, df, 1, 10, 2, 10, std::pow(10, -8));
	std::cout << "result: " << result << std::endl;

	std::cin >> result;
	return 0;
}

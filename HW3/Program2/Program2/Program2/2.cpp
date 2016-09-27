#include <iostream>
#include "Fixed_Point.h"
#include <cmath>

double g(double x) {
	return std::pow(x + 10, 0.25);
}

double g2(double x) {
	return std::pow(x + 10, 0.5) / x;
}

int main(void) {
	double result = fixed_point(g2, 1.0, 10, 0.0000001);
	std::cout << "Result: " << result << std::endl;

	return 0;
}
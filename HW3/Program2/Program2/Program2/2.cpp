// Andrew Sheridan
// Math 5610
// 09/26/15

//Main.cpp
#include <iostream>
#include "Fixed_Point.h"
#include <cmath>

//Test functions found at https://mat.iitm.ac.in/home/sryedida/public_html/caimna/transcendental/iteration%20methods/fixed-point/iteration.html

//First test function.
double g(double x) {
	return std::pow(x + 10, 0.25);
}

//Second test function.
double g2(double x) {
	return std::pow(x + 10, 0.5) / x;
}

int main(void) {
	double result = fixed_point(g, 1.0, 10, std::pow(10, -8));
	std::cout << "Result: " << result << std::endl;
	double result2 = fixed_point(g2, 1.8, 100, std::pow(10, -8));
	std::cout << "Result2: " << result2 << std::endl;

	std::cin >> result;
	return 0;
}
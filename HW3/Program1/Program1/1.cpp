// Andrew Sheridan
// Math 5610
// 09/26/15

#include <iostream>
#include <array>
#include "Bisection.h"

using namespace std;

//First function from example 3.3
double func(double x) { 
	return pow(x, 3) - (30 * x * x) + 2552;
}

//Second function from example 3.3
double func2(double x) {
	return 2.5*sinh(x / 4) - 1;
}

int main() {
	cout << "Testing our bisect method with formulas in Example 3.3." << std::endl;
	double result = bisection(func, 0.0, 20.0, 1 * pow(10, -8));
	cout << "Result 1: " << result << std::endl;
	double result2 = bisection(func2, -10.0, 10.0, 1 * pow(10, -10));
	cout << "Result 2: " << result2 << std::endl;
	cin >> result;

	return 0;
}
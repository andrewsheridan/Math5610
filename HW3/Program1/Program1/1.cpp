#include <iostream>
#include <array>
#include "Bisection.h"

using namespace std;

double func(double x) {
	return pow(x, 3) - (30 * x * x) + 2555;
}

int main() {
	try {
		double result = bisection(func, 0, 20, 1 * pow(10, -8));
		cout << "Result of bisection: " << result << std::endl;
		cin >> result;
	}
	catch (...) {
		cout << "An exception occurred" << std::endl;
	}

	return 0;
}
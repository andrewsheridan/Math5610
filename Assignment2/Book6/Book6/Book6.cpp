#include <iostream>
#include <vector>
#include <cmath>
#include "VectorInput.h"

double doRound(double x, int n) {
	double result = x * pow(10, n);
	result += 0.5;
	result = std::floor(result);
	result *= pow(10, -n);
	return result;
}

std::vector<double> myRound(std::vector<double> list, int n) {
	std::vector<double> result;
	for (int i = 0; i < list.size(); i++) {
		double aThing = doRound(list[i], n);
		result.push_back(aThing);
	}
	return result;
}

int main(void) {
	int n = 0;
	while (true) {
		std::vector<double> list = getVectorInput();
		std::cout << "Please input an integer for your exponent: ";
		std::cin >> n;

		std::vector<double> result = myRound(list, n);
		for (int i = 0; i < result.size(); i++) {
			std::cout << list[i] << ": " << result[i] << std::endl;
		}
	}

	return 0;
}

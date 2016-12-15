#pragma once
#include <cmath>
#include "Complex.h"

double realAbsolute(double number, double approximation) {
	return std::abs(number - approximation);
}

double complexAbsolute(Complex complex, Complex complexApproximation) {
	return (complex - complexApproximation).absoluteError();
}

double realRelative(double number, double approximation) {
	return std::abs(number - approximation) / std::abs(number);
}

double complexRelative(Complex a, Complex b) {
	return (a - b).absoluteError() / a.absoluteError();
}


///Rounds the number x at the exponent n
double round(double x, int n) {
	double result = x * pow(10, n);
	result += 0.5;
	result = std::floor(result);
	result *= pow(10, -n);
	return result;
}

float singlePrecision() {
	float epsilon = 1;
	while (1 + epsilon != 1) {
		epsilon /= 2;
	}
	return epsilon;
}

double doublePrecision() {
	double epsilon = 1;
	while (1 + epsilon != 1) {
		epsilon /= 2;
	}
	return epsilon;
}
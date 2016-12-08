#pragma once
#include <cmath>
#include "Complex.h"

double realAbsolute(double number, double approximation) {
	return std::abs(number - approximation);
}

double imaginaryAbsolute(Complex complex, Complex complexApproximation) {
	return complexAbsolute(complex - complexApproximation);
}

double realRelative(double number, double approximation) {
	return std::abs(number - approximation) / std::abs(number);
}

double imaginaryRelative(Complex a, Complex b) {
	return complexAbsolute(a - b) / complexAbsolute(a);
}


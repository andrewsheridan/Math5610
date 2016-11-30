#pragma once
#include <iostream>
#include <cmath>


double doRound(double x, int n) {
	double result = x * pow(10, n);
	result += 0.5;
	result = std::floor(result);
	result *= pow(10, -n);
	return result;
}
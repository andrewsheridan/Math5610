#pragma once
#include <cmath>

double nested(double data[], unsigned length, double x) {
	double p = data[0];
	for (int i = 1; i < length; i++) {
		p = p*x + data[i];
	}
	return p;
}
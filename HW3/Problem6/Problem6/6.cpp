//main.cpp
#include <iostream>
#include <fstream>
#include "Nested.h"

int main() {
	std::ofstream file;
	file.open("output.txt");
	double list[10] = { 1, -18, 144, -672, 2016, -4032, 5376, -4608, 2304, -512 };
	double step = (2.08 - 1.92) / 161;
	double nestedResults[161];
	unsigned i = 0;
	file << "x \t nested result \t direct result" << std::endl;
	for (double x = 1.92; i < 161; x += step, i++) {
		double nestedResult = nested(list, 10, x);
		nestedResults[i] = nestedResult;

		double directResult = std::pow(x - 2, 9);

		file << x << "\t" << nestedResult << "\t"  << directResult << std::endl;
	}
	file.close();
	return 0;
}
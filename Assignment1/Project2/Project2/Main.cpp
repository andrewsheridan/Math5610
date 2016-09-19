#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>

int main() {
	std::ofstream myfile;
	myfile.open("problem2.txt");
	double stoppingPoint = pow(10, -16);

	std::cout << "h \t" << "approximation \t" << "bigO" << std::endl;
	myfile << "h \t" << "approximation \t" << "bigO" << std::endl;

	for (double h = 1.0; h > stoppingPoint; h *= 0.1) {
		double approximation = (-2 / exp(1)) - (exp(-1 - 2*h) / h) - (1 / (exp(1)*h));
		approximation = abs(approximation);

		double bigO = 2 * h * abs(exp(-1));

		std::cout << h << "\t" << approximation << "\t" << bigO << std::endl;
		myfile << h << "\t" << approximation << "\t" << bigO << std::endl;
	}
	myfile.close();

	myfile.open("problem4.txt");
	stoppingPoint = pow(10, -16);

	std::cout << "h \t" << "approximation \t" << "bigO" << std::endl;
	myfile << "h \t" << "approximation \t" << "bigO" << std::endl;

	for (double h = 1.0; h > stoppingPoint; h *= 0.1) {
		double approximation = cos(1.2) - (sin(1.2+h) - sin(1.2-h))/(2*h);
		double bigO = pow(h, 2) / 6 - cos(1.2);

		std::cout << h << "\t" << approximation << "\t" << bigO << std::endl;
		myfile << h << "\t" << approximation << "\t" << bigO << std::endl;
	}
	myfile.close();

	std::cin >> stoppingPoint;

	return 0;
}
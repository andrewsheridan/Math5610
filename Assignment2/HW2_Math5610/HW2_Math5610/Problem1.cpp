#include "iostream"
#include <stdlib.h>
#include "Complex.h"

void realAbsolute() {
	double number;
	double approximation;
	std::cout << "Input a number: ";
	std::cin >> number;
	std::cout << "Input an approximation: ";
	std::cin >> approximation;

	double result = std::abs(number - approximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void imaginaryAbsolute() {
	double real;
	double imaginary;
	double realApproximation;
	double imaginaryApproximation;
	std::cout << "Input a real number: ";
	std::cin >> real;
	std::cout << "Input an imaginary number: ";
	std::cin >> imaginary;
	std::cout << "Input the real portion of your approximation: ";
	std::cin >> realApproximation;
	std::cout << "Input the imaginary portion of your approximation: ";
	std::cin >> imaginaryApproximation;

	Complex complex(real, imaginary);
	Complex complexApproximation(realApproximation, imaginaryApproximation);

	double result = complexAbsolute(complex - complexApproximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void realRelative() {
	double number;
	double approximation;
	std::cout << "Input a number: ";
	std::cin >> number;
	std::cout << "Input an approximation: ";
	std::cin >> approximation;

	double result = std::abs(number - approximation) / std::abs(number);
	std::cout << "Result: " << result << std::endl << std::endl;
}

Complex imaginaryRelative(Complex a, Complex b) {
	double result = complexAbsolute(a - b) / complexAbsolute(a);
	std::cout << "Result: " << result << std::endl << std::endl;
}

Complex testImaginaryRelative() {
	double real;
	double imaginary;
	double realApproximation;
	double imaginaryApproximation;
	std::cout << "Input a real number: ";
	std::cin >> real;
	std::cout << "Input an imaginary number: ";
	std::cin >> imaginary;
	std::cout << "Input the real portion of your approximation: ";
	std::cin >> realApproximation;
	std::cout << "Input the imaginary portion of your approximation: ";
	std::cin >> imaginaryApproximation;

	Complex complex(real, imaginary);
	Complex complexApproximation(realApproximation, imaginaryApproximation);

	return imaginaryRelative(complex, complexApproximation);
}

int main(void) {
	int input = 1;
	while (input != 0) {
		std::cout << "Please select your desired input algorithm." << std::endl;
		std::cout << "1: Real absolute approximation." << std::endl;
		std::cout << "2: Imaginary absolute approximation." << std::endl;
		std::cout << "3: Real relative approximation." << std::endl;
		std::cout << "4: Imaginary relative approximation." << std::endl;
		std::cout << "0: Exit." << std::endl << std::endl;

		std::cin >> input;

		switch (input) {
		case 1:
			realAbsolute();
			break;
		case 2:
			imaginaryAbsolute();
			break;
		case 3:
			realRelative();
			break;
		case 4:
			imaginaryRelative();
			break;
		default:
			break;
		}
	}

	return 0;
}
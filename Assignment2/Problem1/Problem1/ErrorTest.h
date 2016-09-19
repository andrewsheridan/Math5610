#pragma once

#include <iostream>
#include "Error.h"
#include "Complex.h"

void testRealAbsolute() {
	double number;
	double approximation;
	std::cout << "Input a number: ";
	std::cin >> number;
	std::cout << "Input an approximation: ";
	std::cin >> approximation;

	double result = realAbsolute(number, approximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testRealAbsolute(double number, double approximation) {
	double result = realAbsolute(number, approximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testImaginaryAbsolute() {
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

	double result = imaginaryAbsolute(complex, complexApproximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testImaginaryAbsolute(double real, double imaginary, double realApproximation, double imaginaryApproximation) {
	Complex complex(real, imaginary);
	Complex complexApproximation(realApproximation, imaginaryApproximation);

	double result = imaginaryAbsolute(complex, complexApproximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testRealRelative() {
	double number;
	double approximation;
	std::cout << "Input a number: ";
	std::cin >> number;
	std::cout << "Input an approximation: ";
	std::cin >> approximation;

	double result = realRelative(number, approximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testRealRelative(double number, double approximation) {
	double result = realRelative(number, approximation);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testImaginaryRelative() {
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

	double result = imaginaryRelative(complex, complexApproximation);

	std::cout << "Result: " << result << std::endl << std::endl;
}

double testImaginaryRelative(double real, double imaginary, double realApproximation, double imaginaryApproximation) {
	Complex complex(real, imaginary);
	Complex complexApproximation(realApproximation, imaginaryApproximation);

	return imaginaryRelative(complex, complexApproximation);
}
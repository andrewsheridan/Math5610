#pragma once
#include<cmath>

class Complex
{
public:
	double real;
	double imaginary;

	Complex(double newReal, double newImaginary)
	{
		real = newReal;
		imaginary = newImaginary;
	}

	Complex operator+(Complex c) {
		double realResult = real + c.real;
		double imaginaryResult = imaginary + c.imaginary;
		Complex result(realResult, imaginaryResult);
		return result;
	}

	Complex operator-(Complex c) {
		double realResult = real - c.real;
		double imaginaryResult = imaginary - c.imaginary;
		Complex result(realResult, imaginaryResult);
		return result;
	}

	Complex operator*(Complex c) {
		double realResult = real * c.real - (imaginary * c.imaginary);
		double imaginaryResult = real * c.imaginary + c.real * imaginary;
		Complex result(realResult, imaginaryResult);
		return result;
	}

	Complex operator/(Complex c) {
		double realNumerator = real*c.real - imaginary*c.imaginary;
		double imaginaryNumerator = imaginary*c.real - real*c.imaginary;
		double denominator = (c.real*c.real) - (c.imaginary*c.imaginary);

		return Complex(realNumerator / denominator, imaginaryNumerator / denominator);
	}

	double absoluteError() {
		return std::sqrt(std::pow(real, 2) + std::pow(imaginary, 2));
	}
};
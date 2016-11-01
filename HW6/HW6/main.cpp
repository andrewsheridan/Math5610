//Andrew Sheridan
//Math 5610 
//Written in C++

//main.cpp
#include "Matrix.h";
#include "Vector.h"
#include <iostream>

int main() 
{
	

	//Problem 3: 1-Norm
	int size3 = 4;
	double** matrix3 = CreateMatrix(size3);
	matrix3[size3 - 2][size3 - 1] *= 10;
	double result3 = OneNorm(matrix3, size3);

	std::cout << "Our test matrix and the result of computing the 1-Norm of the matrix." << std::endl;
	PrintMatrix(matrix3, size3);
	std::cout << result3 << std::endl << std::endl;

	//Problem 4: Infinity-Norm
	int size4 = 4;
	double** matrix4 = CreateMatrix(size4);
	matrix4[size4 - 2][size4 - 1] *= 10;
	double result4 = InfinityNorm(matrix4, size4);

	std::cout << "Our test matrix and the result of computing the Infinity-Norm of the matrix." << std::endl;
	PrintMatrix(matrix4, size4);
	std::cout << result4 << std::endl << std::endl;

	int input; 
	std::cin >> input;

	return 0;
}
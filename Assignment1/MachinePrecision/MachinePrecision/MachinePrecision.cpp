// MachinePrecision.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "iostream"
#include "math.h"

int main()
{
	double doubleEpsilon = 1.0; 
	int power = 0;

	while (doubleEpsilon + 1 != 1) {
		power--;
		doubleEpsilon /= 2;
	}

	std::cout << "Epsilon: " << doubleEpsilon << std::endl;
	std::cout << "Power: " << power << std::endl;

	float singleEpsilon = 1.0;
	power = 0;

	while (singleEpsilon + 1 != 1) {
		power--;
		singleEpsilon /= 2;
	}

	std::cout << "Single Epsilon: " << singleEpsilon << std::endl;
	std::cout << "Power: " << power << std::endl;
	std::cin >> singleEpsilon;

	//Problem 2


    return 0;
}


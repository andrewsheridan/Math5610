#include <iostream>
#include "NormTest.h"

int main(void) {
	std::cout << "Testing Euclidean Length" << std::endl;
	std::vector<double> test1;
	test1.push_back(3);
	test1.push_back(4);
	testEuclideanLength(test1);

	std::cout << "Testing Euclidean Length" << std::endl;
	std::vector<double> test2;
	test2.push_back(3.5);
	test2.push_back(4.1);
	test2.push_back(8.4);
	testEuclideanLength(test2);

	std::cout << "Testing manhattan norm" << std::endl;
	std::vector<double> test3;
	test3.push_back(1.1);
	test3.push_back(2.2);
	test3.push_back(-3.3);
	test3.push_back(4.4);
	test3.push_back(-5.5);
	testManhattanNorm(test3);

	std::cout << "Testing Infinity Norm" << std::endl;
	std::vector<double> test4;
	test4.push_back(1.1);
	test4.push_back(2.2);
	test4.push_back(-3.3);
	test4.push_back(4.4);
	test4.push_back(-5.5);
	testInfinityNorm(test4);

	int input = 1;
	while (input != 0) {
		std::cout << "Please select one of the following options: " << std::endl;
		std::cout << "1: Test Euclidean Length" << std::endl;
		std::cout << "2: Test the L1 norm" << std::endl;
		std::cout << "3: Test the infinity norm" << std::endl;
		std::cout << "0: Exit" << std::endl << std::endl;
		std::cin >> input;

		switch (input) {
		case 1:
			testEuclideanLength();
			break;
		case 2:
			testManhattanNorm();
			break;
		case 3:
			testInfinityNorm();
			break;
		case 0:
			break;
		}
	}
	
	return 0;
}
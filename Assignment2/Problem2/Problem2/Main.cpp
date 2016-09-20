#include <iostream>
#include "EuclideanTest.h"

int main(void) {
	int input = 1;
	while (input != 0) {
		std::cout << "Please select one of the following options: " << std::endl;
		std::cout << "1: Test Euclidean Length" << std::endl;
		std::cout << "2: Test the L1 norm" << std::endl;
		std::cout << "3: Test the Linfinity norm" << std::endl;
		std::cout << "0: Exit" << std::endl << std::endl;
		std::cin >> input;

		switch (input) {
		case 1:
			testEuclideanLength();
			break;
		case 2:
			testManhattanNorm();
		case 3:

		case 0:

			break;
		}
		
	}
	
	return 0;
}
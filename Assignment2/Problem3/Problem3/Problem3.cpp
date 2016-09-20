#include<iostream>
#include "DotProduct.h"

int main(void) {
	int input = 1;
	while (input != 0) {
		std::cout << "Please select one of the following options: " << std::endl;
		std::cout << "1: Test Dot Product" << std::endl;
		std::cout << "2: Test Cross Product" << std::endl;
		std::cout << "0: Exit" << std::endl << std::endl;
		std::cin >> input;

		switch (input) {
		case 1:
			testDotProduct();
			break;
		case 2:
			testCrossProduct();
			break;
		default:
			break;
		}
	}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
	return 0;
}
#include<iostream>
#include "DotProduct.h"

int main(void) {
	std::cout << "Testing Dot Product" << std::endl;
	std::vector<double> test1;
	test1.push_back(3);
	test1.push_back(4);
	std::vector<double> test2;
	test2.push_back(10);
	test2.push_back(15);
	testDotProduct(test1, test2);

	std::cout << "Testing Dot Product" << std::endl;
	std::vector<double> test3;
	test3.push_back(3);
	test3.push_back(4);
	test3.push_back(5);
	std::vector<double> test4;
	test4.push_back(0.5);
	test4.push_back(0.4);
	test3.push_back(0.7);
	testDotProduct(test3, test4);

	std::cout << "Testing Cross Product" << std::endl;
	std::vector<double> test5;
	std::vector<double> test6;
	test5.push_back(2);
	test5.push_back(3);
	test5.push_back(4);
	test6.push_back(6);
	test6.push_back(5);
	test6.push_back(1);
	testCrossProduct(test5, test6);

	std::cout << "Testing Cross Product" << std::endl;
	std::vector<double> test7;
	std::vector<double> test8;
	test7.push_back(1.5);
	test7.push_back(2.3);
	test7.push_back(4.1);
	test8.push_back(8.1);
	test8.push_back(2.4);
	test8.push_back(1);
	testCrossProduct(test7, test8);

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
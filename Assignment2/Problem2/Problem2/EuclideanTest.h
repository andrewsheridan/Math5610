#pragma once
#include "Euclidean.h"
#include <iostream>

void testEuclideanLength(std::vector<double> list) {

	std::cout << "Numbers in list: " << std::endl;
	for (int i = 0; i < list.size(); i++) {
		std::cout << list[i] << std::endl;
	}

	double result = euclideanLength(list);
	std::cout << "Result: " << result << std::endl;
}

void testEuclideanLength() {

	std::cout << "Please input some numbers, and input 0 when you're done: " << std::endl;
	double input = 1;
	std::vector<double> list;
	while (input != 0) {
		std::cin >> input;
		list.push_back(input);
	}
	
	double result = euclideanLength(list);
	std::cout << "Result: " << result << std::endl;
}

void testManhattanNorm() {
	std::cout << "Please input some numbers, and input 0 when you're done: " << std::endl;
	double input = 1;
	std::vector<double> list;
	while (input != 0) {
		std::cin >> input;
		list.push_back(input);
	}

	double result = manhattanNorm(list);
	std::cout << "Result: " << result << std::endl;
}
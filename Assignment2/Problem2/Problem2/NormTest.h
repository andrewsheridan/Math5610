#pragma once
#include "Norm.h"
#include <iostream>

std::vector<double> getVectorInput() {
	std::cout << "Please input some numbers, and input 0 when you're done: " << std::endl;
	double input = 1;
	std::vector<double> list;
	while (input != 0) {
		std::cin >> input;
		list.push_back(input);
	}

	return list;
}


void testEuclideanLength(std::vector<double> list) {
	std::cout << "Numbers in list: " << std::endl;
	for (int i = 0; i < list.size(); i++) {
		std::cout << list[i] << std::endl;
	}

	double result = euclideanLength(list);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testEuclideanLength() {
	std::vector<double> list = getVectorInput();
	
	double result = euclideanLength(list);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testManhattanNorm(std::vector<double> list) {
	std::cout << "Numbers in list: " << std::endl;
	for (int i = 0; i < list.size(); i++) {
		std::cout << list[i] << std::endl;
	}

	double result = manhattanNorm(list);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testManhattanNorm() {
	std::vector<double> list = getVectorInput();

	double result = manhattanNorm(list);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testInfinityNorm(std::vector<double> list) {
	std::cout << "Numbers in list: " << std::endl;
	for (int i = 0; i < list.size(); i++) {
		std::cout << list[i] << std::endl;
	}

	double result = infinityNorm(list);
	std::cout << "Result: " << result << std::endl << std::endl;
}

void testInfinityNorm() {
	std::vector<double> list = getVectorInput();

	double result = infinityNorm(list);
	std::cout << "Result: " << result << std::endl << std::endl;
}
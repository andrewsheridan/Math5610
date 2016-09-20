#pragma once
#include<vector>
#include<iostream>

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

std::vector<double> getVectorInput(int size) {
	std::cout << "Please input " << size << " values: " << std::endl;
	double input;
	std::vector<double> list;

	for (int i = 0; i < size; i++) {
		std::cin >> input;
		list.push_back(input);
	}

	return list;
}
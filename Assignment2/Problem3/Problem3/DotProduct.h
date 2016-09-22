#pragma once

#include<cmath>
#include<vector>
#include<iostream>
#include<assert.h>

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

double dotProduct(std::vector<double> one, std::vector<double> two) {
	if (one.size() != two.size()) {
		std::cout << "These vectors are not the same size." << std::endl;
		return -1;
	}

	double sum = 0;

	for (int i = 0; i < one.size(); i++) {
		sum += one[i] * two[i];
	}

	return sum;
}

void testDotProduct() {
	std::vector<double> one = getVectorInput();
	std::vector<double> two = getVectorInput();

	double result = dotProduct(one, two);

	std::cout << "Result: " << result << std::endl;
}

void testDotProduct(std::vector<double> one, std::vector<double> two) {
	std::cout << "Numbers in list one: " << std::endl;
	for (int i = 0; i < one.size(); i++) {
		std::cout << one[i] << std::endl;
	}
	std::cout << "Numbers in list two: " << std::endl;
	for (int i = 0; i < two.size(); i++) {
		std::cout << two[i] << std::endl;
	}
	double result = dotProduct(one, two);
	std::cout << "Result: " << result << std::endl;
}

std::vector<double> crossProduct(std::vector<double> one, std::vector<double> two) {
	std::vector<double> result;
	if (one.size() != two.size()) {
		std::cout << "These vectors are not the same size." << std::endl;
		return result;
	}
	if (one.size() != 3) {
		std::cout << "These vectors are not the correct size." << std::endl;
		return result;
	}

	result.push_back(one[1]*two[2] - one[2]*two[1]);
	result.push_back(one[2] * two[0] - one[0] * two[2]);
	result.push_back(one[0] * two[1] - one[1] * two[0]);
	return result;
}

void testCrossProduct() {
	std::vector<double> one = getVectorInput(3);
	std::vector<double> two = getVectorInput(3);

	std::vector<double> result = crossProduct(one, two);

	std::cout << "Result: " << std::endl;
	std::cout << "<" << result[0] << ", " << result[1] << ", " << result[2] << ">" << std::endl;
}

void testCrossProduct(std::vector<double> one, std::vector<double> two) {
	std::cout << "Numbers in list one: " << std::endl;
	for (int i = 0; i < one.size(); i++) {
		std::cout << one[i] << std::endl;
	}
	std::cout << "Numbers in list two: " << std::endl;
	for (int i = 0; i < two.size(); i++) {
		std::cout << two[i] << std::endl;
	}
	std::vector<double> result = crossProduct(one, two);
	std::cout << "Result: " << std::endl;
	for (int i = 0; i < result.size(); i++) {
		std::cout << result[i] << std::endl;
	}
}

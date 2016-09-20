#pragma once
#include<vector>
#include<cmath>

double euclideanLength(std::vector<double> list) {
	double sum = 0;
	for (int i = 0; i < list.size(); i++) {
		sum += (list[i] * list[i]);
	}
	return std::sqrt(sum);
}

double manhattanNorm(std::vector<double> list) {
	double sum = 0;
	for (int i = 0; i < list.size(); i++) {
		sum += std::abs(list[i]);
	}
	return sum;
}

double infinityNorm(std::vector<double> list) {
	double max = 0;
	for (int i = 0; i < list.size(); i++) {
		if (std::abs(list[i]) > max) {
			max = std::abs(list[i]);
		}
	}
	return max;
}
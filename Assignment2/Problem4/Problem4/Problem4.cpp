#include<iostream>


int main(void) {
	float epsilon = 1;
	while (1 + epsilon != 1) {
		epsilon /= 2;
	}

	std::cout << "Result: " << epsilon << std::endl;
	std::cin >> epsilon;
	return 0;
}
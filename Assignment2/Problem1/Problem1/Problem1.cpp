#include "iostream"
#include <stdlib.h>
#include "Complex.h"
#include "Error.h"
#include "ErrorTest.h"

int main(void) {
	int input = 1;
	while (input != 0) {
		std::cout << "Please select your desired input algorithm." << std::endl;
		std::cout << "1: Real absolute approximation." << std::endl;
		std::cout << "2: Imaginary absolute approximation." << std::endl;
		std::cout << "3: Real relative approximation." << std::endl;
		std::cout << "4: Imaginary relative approximation." << std::endl;
		std::cout << "0: Exit." << std::endl << std::endl;

		std::cin >> input;

		switch (input) {
		case 1:
			testRealAbsolute();
			break;
		case 2:
			testImaginaryAbsolute();
			break;
		case 3:
			testRealRelative();
			break;
		case 4:
			testImaginaryRelative();
			break;
		default:
			break;
		}
	}

	return 0;
}
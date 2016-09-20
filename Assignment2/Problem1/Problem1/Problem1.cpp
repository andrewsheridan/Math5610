#include "iostream"
#include <stdlib.h>
#include "Complex.h"
#include "Error.h"
#include "ErrorTest.h"

int main(void) {
	testRealAbsolute(1.0, 0.99);
	testRealAbsolute(1.0, 1.01);
	testRealAbsolute(100.0, 100.3);
	testRealAbsolute(-1.5, -1.2);
	testRealAbsolute(100.0, 99.99);
	testRealAbsolute(100.0, 99);

	testRealRelative(1.0, 0.99);
	testRealRelative(1.0, 1.01);
	testRealRelative(100.0, 100.3);
	testRealRelative(-1.5, -1.2);
	testRealRelative(100.0, 99.99);
	testRealRelative(100.0, 99);

	testImaginaryAbsolute(3.0, 4.0, 3.0, 4.1);
	testImaginaryAbsolute(3.0, 4.0, 3.3, 4.0);
	testImaginaryAbsolute(3.0, 4.0, 3.2, 4.3);

	testImaginaryRelative(3.0, 4.0, 3.0, 4.1);
	testImaginaryRelative(3.0, 4.0, 3.3, 4.0);
	testImaginaryRelative(3.0, 4.0, 3.2, 4.3);

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
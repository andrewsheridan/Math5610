#pragma once
//Andrew Sheridan
//Math 5610 
//Written in C++

//Vector.h

/// Generates a random vector of size n
// n: The size of the vector
double* CreateVector(unsigned n) {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	double* vector;
	vector = new double[n];
	for (unsigned i = 0; i < n; i++) {
		vector[i] = dis(generator); //Assign each entry in vector to random number between 0 and 1
	}

	return vector;
}

///Creates a vector of all ones of size n
double* CreateOnesVector(unsigned n) {
	double* vector = new double[n];
	for (unsigned i = 0; i < n; i++) {
		vector[i] = 1;
	}
	return vector;
}

///Returns a copy of the input vector
double* CopyVector(double* vector, unsigned size) {
	double* result = new double[size];
	for (unsigned i = 0; i < size; i++) {
		result[i] = vector[i];
	}
	return result;
}

//Finds the entry in V (a size n vector) with the largest magnitude, starting with entry "start".
double FindArrayMax(double* V, unsigned start, unsigned size) {
	double max = 0;
	for (int i = start; i < size; i++) {
		double value = std::abs(V[i]);
		if (value > max) max = value;
	}
	return max;
}

//Finds the index of the value with the largest magnitude in a vector V with size n
int FindMaxIndex(double* V, unsigned n) {
	double max = 0;
	int index = -1;
	for (int i = 0; i < n; i++) {
		double value = std::abs(V[i]);
		if (value > max) {
			max = value;
			index = i;
		}
	}
	return index;
}

///Outputs a size n vector to the console
void PrintVector(double* vector, unsigned size) {
	for (unsigned i = 0; i < size; i++) {
		std::cout << std::setw(13) << std::left << vector[i];
	}
	std::cout << std::endl << std::endl;
}

double DotProduct(double* A, double* B, unsigned int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += A[i] * B[i];
	}
	return sum;
}
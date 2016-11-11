//Andrew Sheridan
//Math 5610 
//Written in C++

#include "Vector.h"
#include <random>
#include <iostream>
#include <iomanip>

Vector::Vector() :size(0){

}

/// Initializes the entries to an empty array of size n
Vector::Vector(unsigned n) : size(n){
	entries = new double[size];
}

///Copy Constructor
Vector::Vector(const Vector &v) : size(v.size) {
	entries = new double[size];
	for (unsigned i = 0; i < size; i++) {
		entries[i] = v.entries[i];
	}
}

Vector::Vector(double* v, unsigned size) {
	entries = new double[size];
	for (unsigned i = 0; i < size; i++) {
		entries[i] = v[i];
	}
}

// Copy Assignment Operator
Vector Vector::operator= (const Vector& v) {
	Vector newVector(v.size);
	for (unsigned i = 0; i < size; i++) {
		newVector[i] = v.entries[i];
	}
	return newVector;
}

//Destructor
Vector::~Vector()
{
}

/// Initializes the entries to values between 0 and 1
void Vector::InitializeRandomEntries() {
	std::mt19937 generator(123); //Random number generator
	std::uniform_real_distribution<double> dis(0.0, 1.0); //Desired distribution

	for (unsigned i = 0; i < size; i++) {
		entries[i] = dis(generator); //Assign each entry in entries to random number between 0 and 1
	}
}

/// Initializes the entries to 1
void Vector::InitializeAllOnes() {
	for (unsigned i = 0; i < size; i++) {
		entries[i] = 1; //Assign each entry in entries to 1
	}
}

//Finds the entry with the largest magnitude, starting with entry "start".
double Vector::FindMaxMagnitude(unsigned start) {
	double max = 0;
	for (unsigned i = start; i < size; i++) {
		double value = std::abs(entries[i]);
		if (value > max) max = value;
	}
	return max;
}

//Finds the index of the value with the largest magnitude
unsigned Vector::FindMaxIndex() {
	double max = 0;
	unsigned index = -1;
	for (unsigned i = 0; i < size; i++) {
		double value = std::abs(entries[i]);
		if (value > max) {
			max = value;
			index = i;
		}
	}
	return index;
}

///Outputs the vector's entries to the console
void Vector::Print() {
	for (unsigned i = 0; i < size; i++) {
		std::cout << std::setw(13) << std::left << entries[i];
	}
	std::cout << std::endl << std::endl;
}

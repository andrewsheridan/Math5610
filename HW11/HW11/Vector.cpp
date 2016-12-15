//Andrew Sheridan
//Math 5610 
//Written in C++
//Vector.cpp

#pragma once
#include "Vector.h"
#include <random>
#include <iostream>
#include <iomanip>

///Default constructor
Vector::Vector() :size(0){ }

/// Initializes the entries to an empty array of size n
Vector::Vector(unsigned n) : size(n){
	entries = new double[size];
	for (unsigned i = 0; i < n; i++) {
		entries[i] = 0;
	}
}

///Copy Constructor
Vector::Vector(const Vector &v) : size(v.size) {
	entries = new double[size];
	for (unsigned i = 0; i < size; i++) {
		entries[i] = v.entries[i];
	}
}

///Copy Constructor
///v: an array of doubles
///n: the size of v
Vector::Vector(double* v, unsigned n) :size(n) {
	entries = new double[size];
	for (unsigned i = 0; i < size; i++) {
		entries[i] = v[i];
	}
}

/// Copy Assignment Operator
Vector Vector::operator= (const Vector& v) {
	size = v.size;
	entries = new double[v.size];
	for (unsigned i = 0; i < v.size; i++) {
		entries[i] = v.entries[i];
	}
	return *this;
}

///Destructor
Vector::~Vector()
{
}

///Inner Product
double operator* (Vector& a, Vector& b) {
	if (a.size != b.size) 
		return NULL;
	double sum = 0;
	for (int i = 0; i < a.size; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

///Multiply each entry by a constant
Vector operator* (Vector& a, double constant) {
	Vector newVector(a.size);
	for (int i = 0; i < a.size; i++) {
		newVector[i] = a[i] * constant;
	}
	return newVector;
}

///Multiply each entry by a constant
Vector operator* (double constant, Vector& a) {
	Vector newVector(a.GetSize());
	for (int i = 0; i < a.GetSize(); i++) {
		newVector[i] = a[i] * constant;
	}
	return newVector;
}

///Subtracts the elements of vector a from the elements of vector b
Vector operator- (Vector& a, Vector& b) {
	if (a.size != b.size) return a;
	Vector newVector(a.size);
	for (int i = 0; i < a.size; i++) {
		newVector[i] = a[i] - b[i];
	}
	return newVector;
}

///Divides the entries of the vector by a constant
Vector operator/ (Vector& a, double constant) {
	Vector newVector(a.GetSize());
	for (int i = 0; i < a.GetSize(); i++) {
		newVector[i] = a[i] / constant;
	}
	return newVector;
}

///Returns a new vector where its entries are the sum of corresponding entries in a and b
Vector operator+ (Vector& a, Vector& b) {
	if (a.size != b.size) {
		std::cout << "Vectors are not same size. Returning first vector." << std::endl;
		return a;
	}
	Vector newVector(a.size);
	for (int i = 0; i < a.size; i++) {
		newVector[i] = a[i] + b[i];
	}
	return newVector;
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

///Finds the entry with the largest magnitude, starting with entry "start".
double Vector::FindMaxMagnitudeStartingAt(unsigned start) {
	double max = 0;
	for (unsigned i = start; i < size; i++) {
		double value = std::abs(entries[i]);
		if (value > max) max = value;
	}
	return max;
}

///Finds the index of the value with the largest magnitude
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

///Computes the L2 norm of the vector
double Vector::L2Norm() {
	double sum = 0;
	for (int i = 0; i < size; i++) {
		sum += (entries[i] * entries[i]);
	}
	return std::sqrt(sum);
}

///Computes the Manhattan norm of the vector
double Vector::ManhattanNorm() {
	double sum = 0;
	for (unsigned i = 0; i < size; i++) {
		sum += std::abs(entries[i]);
	}
	return sum;
}

///Computes the Infinity norm of the vector
double Vector::InfinityNorm() {
	double max = 0;
	for (unsigned i = 0; i < size; i++) {
		if (std::abs(entries[i]) > max) {
			max = std::abs(entries[i]);
		}
	}
	return max;
}

///Outputs the vector's entries to the console
void Vector::Print() {
	for (unsigned i = 0; i < size; i++) {
		std::cout << std::setw(13) << std::left << entries[i];
	}
	std::cout << std::endl << std::endl;
}

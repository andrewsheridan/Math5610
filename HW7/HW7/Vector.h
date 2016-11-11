#pragma once
//Andrew Sheridan
//Math 5610 
//Written in C++

//Vector.h
class Vector {
public:
	const unsigned size; //Size of the vector
	Vector();
	Vector(unsigned n);
	Vector(const Vector &v);
	Vector(double* v, unsigned size);
	Vector operator=(const Vector& v);
	~Vector();

	double& operator[] (unsigned x) { return entries[x]; }

	void InitializeRandomEntries();
	void InitializeAllOnes();

	double FindMaxMagnitude(unsigned start);
	unsigned FindMaxIndex();

	void Print();
	
private:
	double* entries;
};

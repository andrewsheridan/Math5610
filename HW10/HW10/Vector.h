#pragma once
//Andrew Sheridan
//Math 5610 
//Written in C++
//Vector.h

class Vector {
public:
	//Initialization and Destruction
	Vector();
	Vector(unsigned n);
	Vector(const Vector &v);
	Vector(double* v, unsigned size);
	Vector operator=(const Vector& v);
	~Vector();

	void InitializeRandomEntries();
	void InitializeAllOnes();

	//Overloaded Operators
	double& operator[] (unsigned x) { return entries[x]; }
	friend double operator*( Vector& a, Vector& b);
	friend Vector operator*(Vector& a, double constant);
	friend Vector operator*(double constant, Vector& a);
	friend Vector operator+(Vector& a, Vector& b);
	friend Vector operator-(Vector& a, Vector& b);
	friend Vector operator/(Vector& a, double constant);
	friend Vector operator+(Vector& a, Vector& b);

	//Basic Algorithms
	double FindMaxMagnitudeStartingAt(unsigned start);
	unsigned FindMaxIndex();
	double L2Norm();
	double ManhattanNorm();
	double InfinityNorm();

	//Accessing Data
	void Print();
	unsigned GetSize() { return size; }
	void SetSize(unsigned newSize) { size = newSize; }
	
protected:
	unsigned size;

private:
	double* entries; //The stored values of the vector
};

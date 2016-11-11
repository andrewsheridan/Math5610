#pragma once
//Andrew Sheridan
//Math 5610 
//Written in C++

//Vector.h
class Vector {
public:
	Vector();
	Vector(unsigned n);
	Vector(const Vector &v);
	Vector(double* v, unsigned size);
	Vector& operator=(const Vector& v);
	~Vector();

	double& operator[] (unsigned x) { return entries[x]; }
	friend double operator*( Vector& a, Vector& b);
	friend Vector& operator*(Vector& a, double constant);
	friend Vector& operator-(Vector& a, Vector& b);

	void InitializeRandomEntries();
	void InitializeAllOnes();

	double FindMaxMagnitudeStartingAt(unsigned start);
	unsigned FindMaxIndex();

	void Print();
	unsigned GetSize() { return size; }
	void SetSize(unsigned newSize) { size = newSize; }
	
protected:
	unsigned size;

private:
	double* entries;
};

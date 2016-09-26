//#pragma once
//#include<string>
//#include<sstream>
//
//template <typename T, unsigned S>
//unsigned arraysize(const T (&v)[S]) { return S; }
//
//template<typename T>
//struct List {
//	T list[]; 
//
//	List() {
//		list = new T[0];
//	}
//
//	T get(unsigned index) {
//		if (index < arraysize(list)) {
//			return list[index];
//		}
//		else {
//			throw "Index out of range";
//		}
//	}
//
//	void push(T item) {
//		unsigned size = arraysize(list);
//		T newList[] = new T[newsize + 1];
//		for (unsigned i = 0; i < newsize; i++) {
//			newList[i] = list[i];
//		}
//		newList[size] = item;
//		delete[] list;
//		list = newList;
//	}
//
//	std::string toString() {
//		std::stringstream ss;
//		for (unsigned i = 0; i < arraysize(list) - 1; i++) {
//			ss << list[i] << ", ";
//		}
//		ss << list[arraysize(list) - 1];
//		return ss.str();
//	}
//};
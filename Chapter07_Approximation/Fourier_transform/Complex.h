#pragma once
#ifndef COMPLEX_H
#define COMPLEX_H

#include<iostream>

template<typename T>
class Complex {
public:
	Complex(T Re = 0, T Im = 0) : x(Re), y(Im) {}
	~Complex() {}
	Complex(const Complex<T>& C) : x(C.x), y(C.y) {}

	inline T getRe() const { return x; }
	inline T getIm() const { return y; }

	friend std::ostream& operator<<(std::ostream& os, const Complex<T>& C) {
		os << '(' << C.x << ',' << C.y << ")";
		return os;
	}

	Complex<T>& operator=(const Complex<T>& C) { x = C.x; y = C.y; return *this; }

	friend Complex<T> operator+(const Complex<T>& C1, const Complex<T>& C2) {
		return Complex<T>(C1.x + C2.x, C1.y + C2.y);
	}

	friend Complex<T> operator-(const Complex<T>& C1, const Complex<T>& C2) {
		// it can simply like this: return C1 + C2*(-1);
		// but it needs to call operator*, slower. 
		return Complex<T>(C1.x - C2.x, C1.y - C2.y);
	}

	friend Complex<T> operator*(const Complex<T>& C1, const Complex<T>& C2) {
		return Complex<T>(C1.x * C2.x - C1.y * C2.y, C1.x * C2.y + C1.y * C2.x);
	}

	// a*C or C*a				// parameter can exchange order
	friend Complex<T> operator*(const T a, const Complex<T>& C) {
		return Complex<T>(a * C.x, a * C.y);
	}

	// C/a
	Complex<T> operator/(const T a) const {
		if (a == 0) { std::cerr << "Can NOT be divided by zero!" << std::endl; exit(1); }
		return Complex<T>(x / a, y / a);// return (1/a) * (*this);
	}

	// a/C, here parameter cannot exchange order since we also have above member function 
	friend Complex<T> operator/(const T a, Complex<T>& C) {
		if (C.x == 0 && C.y == 0) { std::cerr << "Can NOT be divided by complex zero!" << std::endl; exit(1); }
		return (a / (C.x * C.x + C.y * C.y)) * C.conj();
	}

	// C1/C2
	friend Complex<T> operator/(const Complex<T>& C1, const Complex<T>& C2) {
		if (C2.x == 0 && C2.y == 0) { std::cerr << "Can NOT be divided by complex zero!" << std::endl; exit(1); }
		return C1 * C2.conj() / (C2.x * C2.x + C2.y * C2.y);
	}

	friend double abs(const Complex<T>& C) {
		return sqrt(C.x * C.x + C.y * C.y);
	}

	// return a copy of -C, without changing C itself
	friend Complex<T> operator-(const Complex<T>& C) {
		return Complex<T>(-C.x, -C.y);
	}

	// return a copy of opposite to *this
	Complex<T> oppo() const { return Complex<T>(-x, -y); }

	// render (*this) itself opposite
	Complex<T>& toOppo() { x = -x; y = -y; return *this; }

	// return a copy of its conjugate without changing itself
	Complex<T> conj() const { return Complex<T>(x, -y); }

	// render itself conjugate
	Complex<T>& toConj() { y = -y; return *this; }

private:
	T x;	// real part
	T y;	// imaginary part
};


#endif // !COMPLEX_H

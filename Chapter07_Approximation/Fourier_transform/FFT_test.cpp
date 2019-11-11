#include<cmath>
#include<iostream>


#include"Complex.h"

typedef Complex<double> cp;


int main()
{	
	cp* FFT(const cp* A, const int N);
	cp* IFFT(const cp * x, const int N);

	// N is very large, stack can't hold it, we need to use heap (by new operator)
	const int N = 1024*128;	// 2^17, it will consume 2^17*sizeof(Complex<double>)=2^17*2^4 B=2MB memory
	cp* A = new cp[N];
	for (int i = 0; i < N; ++i)
		A[i] = i;

	cp* x = FFT(A, N);
	std::cout << "After Fourier transform: \n";
	// print 1~4:
	for (int i = 0; i < 4; ++i) {
		std::cout << x[i] << '\n';
	}

	cp* A1 = IFFT(x, N);
	std::cout << "After inverse Fourier transform: \n";
	for (int i = 0; i < 4; ++i) {
		std::cout << A1[i] << '\n';
	}

	delete[] x;
	delete[] A;
	return 0;
}



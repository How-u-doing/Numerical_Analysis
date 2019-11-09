#include<cmath>
#include<complex>
#include<iostream>


#define pi acos(-1)
typedef std::complex<double> cp;


int main()
{
	cp* FFT(cp* A, const int N);
	const int N = 4;
	cp* A = new cp[N];
	for (int i = 1; i < N; ++i)
		A[i] = cp(i,0);	
	cp* x = FFT(A, N);

	std::cout << "After Fourier transform: \n";
	for (int i = 0; i < N; ++i) {
		std::cout << x[i] << '\n';
	}

	std::cout << x[1]+x[2] << '\n';
	delete[] x;
	delete[] A;
	return 0;
}


// Fast Fourier transform
cp* FFT(cp* A, const int N)
{
	int p = int(log2(N));
	if (p != log2(N)) {
		std::cerr << "log2(N) must be an integer!\n";
		exit(1);	
	}

	cp* A0 = new cp[N];
	// duplicate A to evade to change it
	for (int i = 0; i < N; ++i)
		* (A0 + i) = *(A + i);
	
	cp* w = new cp[N / 2 - 1];
	for (int j = 0; j <= N / 2 - 1; ++j)
		w[j] = cp(cos(-2 * pi * j / N), sin(-2 * pi * j / N));

	cp* A1 = new cp[N];
	for (int q = 1; q <= p; ++q) {
		if (q % 2 == 1) {
			for (int k = 0; k < pow(2, p - q); ++k)
				for (int j = 0; j < pow(2, q - 1); ++j) {
					A1[k * int(pow(2, q)) + j] = A0[k * int(pow(2, q - 1)) + j] + A0[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))];
					A1[k * int(pow(2, q)) + j + int(pow(2, q - 1))] = (A0[k * int(pow(2, q - 1)) + j] -
						A0[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))]) * w[k * int(pow(2, q - 1))];
				}
		}
		else
			for (int k = 0; k < pow(2, p - q); ++k)
				for (int j = 0; j < pow(2, q - 1); ++j) {
					A0[k * int(pow(2, q)) + j] = A1[k * int(pow(2, q - 1)) + j] + A1[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))];
					A0[k * int(pow(2, q)) + j + int(pow(2, q - 1))] = (A1[k * int(pow(2, q - 1)) + j] -
						A1[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))]) * w[k * int(pow(2, q - 1))];
				}
	}
	
	// It's very weird that it throws an exception (Stack.exe has triggered a breakpoint.)
	// when I attempt to release w. ??
	delete[] w;

	if (p % 2 == 0) {
		delete[] A1;
		return A0;
	}		
	else {
		delete[] A0;
		return A1;			
	}
}

#include<cmath>
#include"Complex.h"

#define pi acos(-1)

typedef Complex<double> cp;


// Fast Fourier Transform
cp* FFT(const cp* A, const int N)
{
	int p = int(log2(N));
	if (p != log2(N)) {
		std::cerr << "log2(N) must be an integer!\n";
		exit(1);
	}

	cp* A0 = new cp[N];
	// duplicate A
	for (int i = 0; i < N; ++i)
		*(A0 + i) = *(A + i);

	cp* w = new cp[N / 2];
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


// Inverse Fast Fourier Transform
cp* IFFT(const cp* x, const int N)
{
	int p = int(log2(N));
	if (p != log2(N)) {
		std::cerr << "log2(N) must be an integer!\n";
		exit(1);
	}

	cp* x0 = new cp[N];
	// x0=x/N
	for (int i = 0; i < N; ++i)
		*(x0 + i) = *(x + i) / N;

	cp* w = new cp[N / 2];
	for (int j = 0; j <= N / 2 - 1; ++j)
		w[j] = cp(cos(2 * pi * j / N), sin(2 * pi * j / N));

	cp* x1 = new cp[N];
	for (int q = 1; q <= p; ++q) {
		if (q % 2 == 1) {
			for (int k = 0; k < pow(2, p - q); ++k)
				for (int j = 0; j < pow(2, q - 1); ++j) {
					x1[k * int(pow(2, q)) + j] = x0[k * int(pow(2, q - 1)) + j] + x0[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))];
					x1[k * int(pow(2, q)) + j + int(pow(2, q - 1))] = (x0[k * int(pow(2, q - 1)) + j] -
						x0[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))]) * w[k * int(pow(2, q - 1))];
				}
		}
		else
			for (int k = 0; k < pow(2, p - q); ++k)
				for (int j = 0; j < pow(2, q - 1); ++j) {
					x0[k * int(pow(2, q)) + j] = x1[k * int(pow(2, q - 1)) + j] + x1[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))];
					x0[k * int(pow(2, q)) + j + int(pow(2, q - 1))] = (x1[k * int(pow(2, q - 1)) + j] -
						x1[k * int(pow(2, q - 1)) + j + int(pow(2, p - 1))]) * w[k * int(pow(2, q - 1))];
				}
	}

	delete[] w;

	if (p % 2 == 0) {
		delete[] x1;
		return x0;
	}
	else {
		delete[] x0;
		return x1;
	}
}

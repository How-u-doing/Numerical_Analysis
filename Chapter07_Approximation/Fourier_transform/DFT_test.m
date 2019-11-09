% Discrete Fourier Transform test
A=1:2^10;
x=DFT(A);
A1=IDFT(x);

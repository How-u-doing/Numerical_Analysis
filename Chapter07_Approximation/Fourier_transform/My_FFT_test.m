% FFT and IFFT test
A=1:2^15;   % in such no. of input data, use DFT will be out of 
            % memory (it needs 2^30 times complex multiplication)
% Try complex data: A=1i*(1:2^3);  

N=length(A);
x=My_FFT(A);
% compare with system provided fft
% y=fft(A);
% z=x-y;
% z1=z(1:10).'


% compare A with A1
A1=My_IFFT(x);
% see part of A1 randomly
A2=A1(80:90).'


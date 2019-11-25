% Romberg integration test
f=@(x) sin(x);
a=0;
b=pi;
[T,k]=Romberg(f,a,b)
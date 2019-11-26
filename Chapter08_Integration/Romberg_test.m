% Romberg integration test
f=@(x) sin(x);
a=0;
b=pi;
format long
% for maxIt=1:6
%     fprintf('maxIt=%d\n',maxIt)
%     [T,k]=Romberg(f,a,b,maxIt)
% end
[T,k]=Romberg(f,a,b)
[S,k]=Romberg_Simpson(f,a,b)


f=@(x) sin(x);
a=0;
b=pi;
n=2^10;
format long;
method="Simpson";
T=CompositeInt(f,a,b,n,method)
method="trapezoid";
T=CompositeInt(f,a,b,n,method)
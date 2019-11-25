f=@(x) sin(x);
a=0;
b=pi;
n=4;
%method="trapezium";
method="Simpson";

T=CompositeInt(f,a,b,n,method)
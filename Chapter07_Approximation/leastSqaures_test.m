% Least Squares Approximation test
% Compare & play around with APPS->Curve Fitting
x=linspace(1,2,10);
N=length(x);
y=exp(x);         % original data
y1=y+rand(1,N);   % add interference 
n=3;              % approximation degree
a=leastSqaures(x,y1,n);
p=flip(a);
yls=polyval(p,x);
plot(x,y1,'o',x,yls,'r--');
xlabel('x');ylabel('y');
legend('Captured Data','Least Squares')



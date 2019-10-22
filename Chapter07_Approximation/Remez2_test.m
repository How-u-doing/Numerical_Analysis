% Remez's 2nd algorithm test
% 10170437 Mark Taylor
f=@(x)log(x);
a=0.5;b=5;
% Try: --> play around with these input arguments f,a,b,x

n=3; % <-- change degree here
x=linspace(a,b,n+2); % evenly spaced points
%x=[1.1 2.2 3.3 4.4 5]; % <-- or you provide some random test points (sorted) here

y=f(x);
[p,Enf,x,k]=Remez2(f,x,y,a,b)
xx=a:.01:b;
ypn=polyval(p,xx);figure
plot(xx,f(xx),'g',xx,ypn,'--r')
legend('Original Function','Optimal Polynomial')

% If you see the error function having at least n+2 alternative signs of 
% max values (they are equal in absolute values), it means this is the 
% very optimal function we are looking for.
max_err=polyval(p,x)-f(x)  % campare to the summit of the graph to verify
yerr=polyval(p,xx)-f(xx); figure;
plot(xx,yerr)
legend('Error function')



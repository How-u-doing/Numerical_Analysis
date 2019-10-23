% Remez's 2nd algorithm test
% 10170437 Mark Taylor
f=@(x)sin(1+x.^2);
% f=@(x)x.*sin(1+x.^2);
% f=@(x)x.^5; % when n>=5, error<=1e-14
a=-2;b=2;
% Try: --> play around with these input arguments f,a,b,x

n=4; % <-- change degree here
x=linspace(a,b,n+2); % evenly spaced points
%x=[1.1 2.2 3.3 4.4 5]; % <-- or you provide some random test points (sorted) here

y=f(x);
[p,Enf,x,k]=Remez2(f,x,y,a,b)
xx=a:.01:b;
ypn=polyval(p,xx);figure
plot(xx,f(xx),'g',xx,ypn,'--r')
legend('Original Function','Optimal Polynomial')

% If you see the error function having at least n+2 points with alternative 
% signs of max values (equal in absolute values) or within a tiny order of
% magnitude, for example 1e-8 or smaller (? I think 1e-4 may already be cool)
% it means this is the vey optimal polynomial we are looking for.
max_err=polyval(p,x)-f(x)  % campare to the summit of the graph to verify
yerr=polyval(p,xx)-f(xx); figure;
plot(xx,yerr)
legend('Error function')



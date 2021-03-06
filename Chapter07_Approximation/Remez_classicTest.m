% Remez's 1st algorithm test on 1/(1+x^2)
% 10170437 Mark Taylor
f=@(x) 1./(1+x.^2);
a=-2;b=2;
% Try: --> play around with these input arguments f,a,b,x

n=7; % <-- change degree here
x=linspace(a,b,n+2); % evenly spaced points
%x=[-2,-1,-0.8,-0.5,0,0.3,0.8,1.3,1.7]; % <-- or you provide some random test points (sorted) here

syms t;y=f(x);
xx=a:.01:b;
x0=linspace(a,b,n+1); % given points for Newton interpolating polynomial of degree n
y0=f(x0);
ynew=NI(x0,y0,xx);

[pn,Enf,x,k]=Remez(f,x,y,a,b)
ypn=subs(pn,t,xx);figure
plot(xx,f(xx),'g',xx,ypn,'--r',xx,ynew,'--b')
legend('Original Function','Optimal Polynomial','Newton''s P7')
title('7th-degree Optimal Polynomial')

% If you see the error function having at least n+2 alternative signs of 
% max values (they are equal in absolute values), it means this is the 
% very optimal function we are looking for.
err=pn-f(t);
max_err=vpa(subs(err,t,x),5)  % campare to the summit of the graph to verify
yerr=subs(err,t,xx); figure;
plot(xx,yerr)
legend('Error function')


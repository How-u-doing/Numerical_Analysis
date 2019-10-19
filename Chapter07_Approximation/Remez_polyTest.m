% Remez's 1st algorithm's test on polynomials over interval [a,b].
% Of course, when the no. of test points is no less than the degree of f(x) 
% plus 2, we can get goddamn perfect approximation: f(x) itself. But we can 
% use other easier & more effective algorithms to determine such polynomial.
% Our interest is, however, to try to find an optimal polyniomial of lower
% degree that can best approximate it. 
% Let f(x)=an*x^n+a(n-1)*x^(n-1)+...+a1*x+a0, an~=0, then its (n-1)th degree
% optimal polynomial is: g(x)=an*q(x)+a(n-1)*x^(n-1)+...+a1*x+a0, where q(x)
% is the (n-1)th degree optimal polynomial of x^n. And another useful and
% important conclusion is: En(f)=|an|*((b-a)/2)^n * 2^(1-n). Check it!
% 10170437 Mark Taylor
f=@(x)3*x.^5-5*x.^2+9; % En(f)=3*((3+1)/2)^5 * 2^(-4)=6
% f=@(x)x.^5; % get our q(x) and *3, then -5*x.^2+9 to campare above result
a=-1;b=3;
% Try: --> play around with these input arguments f,a,b,x

n=4; % <-- change degree here
x=linspace(a,b,n+2); % evenly spaced test points
%x=[1.1 2.2 3.3 4.4 5]; % <-- or you provide some random test points (sorted) here

syms t;y=f(x);
[pn,Enf,x,k]=Remez(f,x,y,a,b)
xx=a:.01:b;
ypn=subs(pn,t,xx);figure
plot(xx,f(xx),'g',xx,ypn,'--r')
legend('Original Function','Optimal Polynomial')

% If you see the error function having at least n+2 alternative signs of 
% max values (they are equal in absolute values), it means this is the 
% very optimal function we are looking for.
err=pn-f(t);
max_err=vpa(subs(err,t,x),5)  % campare to the summit of the graph to verify
yerr=subs(err,t,xx); figure;
plot(xx,yerr)
legend('Error function')


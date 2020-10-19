% ===============================================================
% Solve a linear two-point boundary-value problem
%                   -(p(x)u'(x))'+q(x)u(x)=f(x),    a<=x<=b,
% with initial conditions
%                           u(a)=0, u'(b)=0.
% ===============================================================

% ***************************************************************
% Test 1
% f=@(x)pi^2/4*sin(pi/2*x);
% u=@(x)sin(pi/2*x);

% Test 2
f=@(x)pi^2*sin(pi*x);
u=@(x)sin(pi*x)+pi*x;
% ***************************************************************

p=@(x)x.^0; % p(x)=1
q=@(x)x.*0; % q(x)=0
a=0;
b=3;
% b=1;
N=[3,5,10];

% Since the dimensions are different, we can't
% simply use a matrix for storage. List?
uc1=PQRR(f,p,q,a,b,N(1)); 
uc2=PQRR(f,p,q,a,b,N(2)); 
uc3=PQRR(f,p,q,a,b,N(3)); 

% output coef vector of approx solution in command window
disp(uc3) 

% plot graphs of exact solution & approximate solutions
tt=linspace(a,b,100+1);

% approx solutions
u1=PQRR_intpol(uc1,tt,a,b);
u2=PQRR_intpol(uc2,tt,a,b);
u3=PQRR_intpol(uc3,tt,a,b);

figure
plot(tt,u(tt),'--b',tt,u1,'r',...
    tt,u2,'g',tt,u3,'m')
legend('exact',sprintf('N=%d',N(1)),sprintf('N=%d',N(2))...
    ,sprintf('N=%d',N(3)),'Location','northwest')



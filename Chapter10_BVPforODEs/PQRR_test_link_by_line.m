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
t1=linspace(a,b,2*N(1)+1);
t2=linspace(a,b,2*N(2)+1);
t3=linspace(a,b,2*N(3)+1);
tt=linspace(a,b,100+1);

% approx solutions
u1=[0;uc1]; % use line instead, NO interpolation
u2=[0;uc2];
u3=[0;uc3];

figure
plot(tt,u(tt),'--b',t1,u1,'r',...
    t2,u2,'g',t3,u3,'m')
legend('exact',sprintf('N=%d',N(1)),sprintf('N=%d',N(2))...
    ,sprintf('N=%d',N(3)),'Location','northwest')



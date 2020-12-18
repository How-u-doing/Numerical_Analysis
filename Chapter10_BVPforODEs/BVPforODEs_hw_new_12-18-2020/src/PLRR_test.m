% ===============================================================
% Solve a linear two-point boundary-value problem
%                   -(p(x)u'(x))'+q(x)u(x)=f(x),    a<=x<=b,
% with initial conditions
%                           u(a)=0, u'(b)=0.
% ===============================================================

% ***************************************************************
% Test 1
isTest1 = true;
f=@(x)pi^2/4*sin(pi/2*x);
u=@(x)sin(pi/2*x);

% Test 2
% isTest1 = false;
% f=@(x)pi^2*sin(pi*x);
% u=@(x)sin(pi*x)+pi*x;
% ***************************************************************

p=@(x)x.^0; % p(x)=1
q=@(x)x.*0; % q(x)=0
a=0;
b=3;
% b=1;
N=[3,5,10]; 
% N=2*[3,5,10]; % compare with quadratic element

% Since the dimensions are different, we can't
% simply use a matrix for storage. List?
uc1=PLRR(f,p,q,a,b,N(1)); 
uc2=PLRR(f,p,q,a,b,N(2)); 
uc3=PLRR(f,p,q,a,b,N(3)); 

% output coef vector of approx solution in command window
disp(uc3) 

% plot graphs of exact solution & approximate solutions
t1=linspace(a,b,N(1)+1);
t2=linspace(a,b,N(2)+1);
t3=linspace(a,b,N(3)+1);
tt=linspace(a,b,100+1);

% approx solutions
u1=PLRR_intpol(t1,[0;uc1],tt); % Note that these interpolating points are also on
u2=PLRR_intpol(t2,[0;uc2],tt); % those line segments, which means we can actually
u3=PLRR_intpol(t3,[0;uc3],tt); % skip this phase, and replace with [0;uc] instead.

figure
plot(tt,u(tt),'--b',tt,u1,'r',...
    tt,u2,'g',tt,u3,'m')
legend('exact',sprintf('N=%d',N(1)),sprintf('N=%d',N(2))...
    ,sprintf('N=%d',N(3)),'Location','northwest')

if isTest1 && b==3
    ylim([-1.5,1.5])
end




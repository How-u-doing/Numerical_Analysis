% Remez's 2nd algorithm£¬to find an optimal polynomial that can best approx f(x)
% Note that this function has to calculate & compare massive function values
% to determine the max error points, which makes k very large, whereas the
% measure adopted in Remez's 1st algorithm (use syms and solve differential
% equations) can elude it. But it also costs a lot for the latter to solve
% differential equations!
% 10170437 Mark Taylor
function [p,Enf,x,k]=Remez2(f,x,y,a,b,epsilon,maxIt)
%INPUT:
    % f: the given function supporting vector input, e.g. f=@(x)sin(x)
    % x: a vecter containing n+2 initial test points
    % y: corresponding function values of x, i.e. f(x)
    % [a,b]: the given interval
    % epsilon: the given tolerance, stopping criterion
    % maxIt: maximum iteration
%OUTPUT:
    % p: coefficients of the optimal polynomial
    % Enf: maximum error
    % x: determined alternative points
    % k: no. of iterations
    
% set default input arguments
if nargin<7
    maxIt=100;
    if nargin<6
        epsilon=1e-6;
        if nargin<5
            error('Error! Insufficient input arguments!\n')
        end
    end
end
if isrow(x)
    x=x.';
end
if isrow(y)
    y=y.';
end
n=length(x)-2; % the optimal polynomial to f(x) of degree at most n
A=zeros(n+2);  % coefficient matrix
A(:,2)=ones(n+2,1);     % coefficient of a0 (=1)
A(1,1)=1;
for i=2:(n+2)    
    A(i,1)=-1*A(i-1,1); % alternative sign for sigma (En(f) or -En(f))  
    % or A(i,1)=1-2*i%2 (mod(i,2)) , i=1,2,...,n+2
end
for j=1:n
    A(:,j+2)=x.^j;    
end
% X0=[sigma, a0, a1, ... , an].' 
X0=A\y;        % initial solutions

p=flip(X0(2:end).');   % coefficients of Pn(x)
eta=zeros(1,n+3);      % n+1 zeros of error function plus a & b
eta(1)=a; eta(n+3)=b;
for k=1:maxIt
    % find n+1 zeros of error function Pn(x)-f(x) on the interval [a,b]
    for i=2:n+2
        eta(i)=FalsePosition(f,p,x(i-1),x(i));
    end
    
    % Let's adopt a seemingly dumb method to find the max error point
    for i=1:n+2
        % this job costs lots of time
        dx=eta(i):1e-4:eta(i+1);
        [M,idx]=max(abs(f(dx)-polyval(p,dx)));
        x(i)=dx(idx);   % update x
    end
    
    % update A, y, X0 & p
    for j=1:n
        A(:,j+2)=x.^j;  
    end 
    y=f(x);
    X1=A\y;
    p=flip(X1(2:end).');
    
    % tell when to stop
    stop=true;
    for i=0:n
        if abs(X1(i+2)- X0(i+2))>epsilon
            stop=false;
            break;
        end
    end
    if stop==true
        Enf=X1(1);
        return;
    else % update X0
        X0=X1;
    end
end
fprintf('Maximum iteration exceeded£¡The results in the last iteration as followed:\n')
Enf=X1(1);
end


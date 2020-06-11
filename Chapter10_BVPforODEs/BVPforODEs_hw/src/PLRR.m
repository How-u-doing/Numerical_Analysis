% Piecewise Linear Rayleigh-Ritz method
% ===============================================================
% Solve a linear two-point boundary-value problem
%                   -(p(x)u'(x))'+q(x)u(x)=f(x),    a<=x<=b,
% with initial conditions
%                           u(a)=0, u'(b)=0.
% ===============================================================
% See the theory part at <reference.mlx>
function uc=PLRR(f,p,q,a,b,N)
% INPUT:
%   f: f(x)
%   p: p(x); optional, by default p(x)=1
%   q: q(x); optional, by default q(x)=0
%   a, b: interval [a,b]; optional, by default [a,b]=[0,1]
%   N: # of evenly spaced intervals; optional, by default N=10
% OUTPUT:
%   uc: coef vector of approx solution

% ===============================================================
% *****  NOTE that f, p, q must be element-wise functions  *****
% e.g. f(x)=@(x)x.^2+1 (NOT x^2+1)
% ===============================================================

% set default args
if nargin<6
    N=10;
    if nargin<5
        b=1;
        if nargin<4
            a=0;
            if nargin<3
                q=@(x)x.*0;
                if nargin<2
                    p=@(x)x.^0;
                    if nargin<1                 
                        error('Error! Function f is not set.')
                    end
                end
            end
        end
    end
end
if b<a
    tmp=a; a=b; b=tmp;
end

h=(b-a)/N;  % step length. no check for N, hhh
A=zeros(N); % coef matrix, symmetric tridiagonal
c=zeros(N,1); % right-hand side (RHS) vector
xx=linspace(a,b,N+1);

% determine tridiagonal elements of A, and right-side vector c
for i=1:N-1
    % Or use MATLAB func <integral> instead of <Gaussquad>.
    % --> <Ctrl + F> --> Replace All. Then go to PLRR_test
    % script file & run it again, it shall get same graph.
    A(i,i)=Gaussquad(p,xx(i),xx(i+2))+...
        Gaussquad(@(x)(x-xx(i)).^2.*q(x),xx(i),xx(i+1))+...
        Gaussquad(@(x)(xx(i+2)-x).^2.*q(x),xx(i),xx(i+1));    
    A(i,i+1)=-Gaussquad(p,xx(i+1),xx(i+2))+...
        Gaussquad(@(x)(x-xx(i+1)).*(xx(i+2)-x).*q(x),xx(i+1),xx(i+2));
    
    c(i)=Gaussquad(@(x)(x-xx(i)).*f(x),xx(i),xx(i+1))+...
        Gaussquad(@(x)(xx(i+2)-x).*f(x),xx(i+1),xx(i+2));   
end
A(N,N)=Gaussquad(p,xx(N),xx(N+1))+...
       Gaussquad(@(x)(x-xx(N)).^2.*q(x),xx(N),xx(N+1));
c(N)=Gaussquad(@(x)(x-xx(N)).*f(x),xx(N),xx(N+1));

A=A+diag(diag(A,1),-1); % let A_{i,i-1}=A_{i-1,i}
A=A/h^2;
c=c/h;

% uc=A\c; % coef vector of approx solution
% or use user-defined function <GauEli>
uc=GauEli(A,c);

end




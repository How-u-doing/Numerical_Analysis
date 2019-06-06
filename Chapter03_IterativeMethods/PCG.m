% Preconditioned Conjugate Gradient Method
% 10170437 Mark Taylor

function [x, r, k]=PCG(A, b, CI, tol, N, x_0)
% INPUT:
%   A: coefficient matrix which is symmetric & positive definite,
%   b: right hand side vector, 
%   CI: the inverse of matrix C
%   tol: tolerance, 
%   N: maximum number of iterations,
%   x_0: initial approximation(by default,x_0=zeros(n,1)).
% OUTPUT:
%   x: approximation solution vector,
%   r: residual vector,   
%   k: the number of iterations.

% The complete procedure to check if A is symmetric & positive definite is: 
% A~=A.'  &&  allPositive(eig(A))==true. For executive efficiency, however, 
% we just check whether it is symmeric.
if  isequal(A,A.')==false
    error('A must be symmetric & positive definite!')
end

n=size(A,2);
% Set default initializations
if nargin<6
    x_0=zeros(n,1);
    if nargin<5
        N=1000;
        if nargin<4
            tol=1e-6;
            if nargin<3
                error('Too few input arguments!')
            end
        end
    end 
end

% Initialization
x=x_0;
r=b-A*x_0;
w=CI*r;
p=CI.'*w;
t=w.'*w;

k=0;

while k<=N
                            
    if norm(r,inf)<tol || norm(p,inf)<tol
        return;
    end
    
    u=A*p;
    alpha=t/(p.'*u);
    x=x+alpha*p;
    r=r-alpha*u;
    w=CI*r;
    s=w.'*w;
    
    beta=s/t; 
    p=CI.'*w+beta*p;
    t=s;
    
    k=k+1;

end

% The number of iterations was exceeded.
k=k-1;% k=N
fprintf('\nCannot compute the approximate solution vector x within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('The last iterative approximate solution vector x is as followed:\n');

end






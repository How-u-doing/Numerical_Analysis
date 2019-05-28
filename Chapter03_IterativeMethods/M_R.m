% Minimum Residual method
% 10170437 Mark Taylor

function [x, k]=M_R(A, b, tol, N, x_0)
% INPUT:
%   A: coefficient matrix which is symmetric & positive definite,
%   b: right hand side vector, 
%   tol: tolerance, 
%   N: maximum number of iterations,
%   x_0: initial approximation(by default,x_0=zeros(n,1)).
% OUTPUT:
%   x: approximation solution vector,   
%   k: the number of iterations.

% The complete procedure to check if A is symmetric & positive definite is: 
% A~=A.'  &&  allPositive(eig(A))==true. For executive efficiency, however, 
% we just check whether it is symmeric.
if  isequal(A,A.')==false
    error('A must be symmetric & positive definite!')
end

n=size(A,2);
if nargin<5
    x_0=zeros(n,1); % set default initial approximation
end

x=x_0;
r=b-A*x_0;
k=1;
while k<=N
    beta=(A*r).'*r/((A*r).'*(A*r));
    x=x+beta*r;
    r=b-A*x;% <-> r=r-beta*A*r; 
    
    if norm(r,inf)<tol
        return;
    end
    
    k=k+1;

end

% The number of iterations was exceeded.
k=k-1;% k=N
fprintf('\nCannot compute the approximate solution vector x within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('The last iterative approximate solution vector x is as followed:');

end






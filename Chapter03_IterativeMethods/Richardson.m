% Richardson's Relaxation Method 
% 10170437 Mark Taylor
function [x, k]=Richardson(A, b, w, tol, N, X0)
% INPUT:
%   A: coefficient matrix,
%   b: right hand side vector, 
%   w: parameter omega,
%   tol: tolerance, 
%   N: maximum number of iterations,
%   X0: initial approximation(by default,X0=zeros(n,1)).
% OUTPUT:
%   x: approximation solution vector,   
%   k: the number of iterations.

[m,n]=size(A);
if m~=n
    error('A must be square!')
end

% Set default initializations
if nargin<6
    X0=zeros(n,1);
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


k=1;
x=X0;
while k<=N
    
    x=X0+w*(b-A*X0);%Here X0 can be replaced by x        
    
    if norm(x-X0,inf)<tol
        return;
    end
    
    X0=x;
    k=k+1;

end

% The number of iterations was exceeded.
k=k-1;% k=N
fprintf('\nCannot compute the approximate solution vector x within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('The last iterative approximate solution vector x is as followed:\n');

end
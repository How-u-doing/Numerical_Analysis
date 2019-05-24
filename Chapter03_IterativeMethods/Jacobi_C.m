% Jacobi Iteration to solve linear systems of algebra equations Ax=b.
% x_(k)=inv(D)(L+U)x_(k-1)+inv(D)b, 
% or x_(k)=(I-inv(D)A)x_(k-1)+inv(D)b, k=1,2,3, ...
% where D is a diagonal matrix whose diagonal entries are those of A,
% -L is the strictly lower-triangular part of A, and -U is the strictly 
% upper-triangular part of A. (A=D-L-U)
% 10170437 Mark Taylor

function [x, k]=Jacobi_C(A, b, tol, N, X0)
% Component-wise manner of Jacobi's Method

% INPUT:
%   A: coefficient matrix
%   b: right hand side vector 
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

% Rearrange A in terms of row so as to make diagnonal entries' absolute 
% value greater than zero and as large as possible to speed convergence.
[isSuccessful, isReordered, Aug]=reordering([A, b]);
A=Aug(:,1:n);
b=Aug(:,n+1);

if isSuccessful==true
    if isReordered==true
        disp('A reordering of the equations is performed successfully!');
        disp('The reordered A is as followed:');
        A
%   else
%       No need to reorder A.  
    end
else
    error('Reordering failure to make all diagonal entries greater than 1.0e-6!')
end

if nargin<5
    X0=zeros(n,1); % set default initial approximation
end

% #######  This part is not necessary in practical computation.
inv_D=diag(1./diag(A));
% Though it is never a cost-effecient way to tell whether iteration matrix 
% converges via using eig function which consumes relatively large amount
% of memory & time. Temporarily I can't figure out a better approach.
if max(abs(eig(eye(n)-inv_D*A)))>=1 
    error('Iterative matrix does not converge!')
end
% #######

k=1;
x=X0;
while k<=N
    for i=1:n
        x(i)=-1/A(i,i)*(A(i,1:i-1)*X0(1:i-1)+A(i,i+1:n)*X0(i+1:n)-b(i));
%                                  ^
%                                  | 
% Jacobi's Method uses these previous values 
        
    end
    
    if norm(x-X0,inf)<tol
%   or norm(x-X0,inf)/norm(x,inf)<tol
        return;
    end
    
    X0=x;
    k=k+1;

end

% The number of iterations was exceeded.
k=k-1;% k=N
fprintf('\nCannot compute the approximate solution vector x within %d iterations in the tolerance of %d!\n',N,tol);
fprintf('The last iterative approximate solution vector x is as followed:');

end






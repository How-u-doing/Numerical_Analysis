% Gauss-Seidel Iteration to solve linear systems of algebra equations Ax=b.
% x_(k)=inv(D)(Lx_(k)+Ux_(k-1)+b), 
% or x_(k)=inv(D-L)Ux_(k-1)+inv(D-L)b, k=1,2,3, ...
% where D is a diagonal matrix whose diagonal entries are those of A,
% -L is the strictly lower-triangular part of A, and -U is the strictly 
% upper-triangular part of A. (A=D-L-U)
% 10170437 Mark Taylor

function [x, k]=GaussSeidel_C(A, b, tol, N, x_0)
% Component-wise form of Gauss-Seidel's Method

% INPUT:
%   A: coefficient matrix
%   b: right hand side vector 
%   tol: tolerance, 
%   N: maximum number of iterations,
%   x_0: initial approximation(by default,x_0=zeros(n,1)).
% OUTPUT:
%   x: approximation solution vector,   
%   k: number of iterations.

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
    x_0=zeros(n,1); % set default initial approximation
end

inv_D=diag(1./diag(A));

% Though it is never a cost-effecient way to tell whether iteration matrix 
% converges via using eig function which consumes relatively large amount
% of memory & time. Temporarily I can't figure out a better approach.
if max(abs(eig(eye(n)-inv_D*A)))>=1 
    error('Iteration matrix does not converge!')
end

k=1;
XP=x_0; % previous iteration of x
x=zeros(n,1);
while k<=N
    for i=1:n
        x(i)=-1/A(i,i)*(A(i,1:i-1)*x(1:i-1)+A(i,i+1:n)*XP(i+1:n)-b(i));
%                                  ^
%                                  | 
% Gauss-Seidel's Method uses these most recently calculated values 
    end
    
    if norm(x-XP,inf)/norm(x-x_0,inf)<tol
        return;
    end
    
    XP=x;
    k=k+1;

end

error('Cannot compute approximation solution vector x within %d iterations!',N)

end






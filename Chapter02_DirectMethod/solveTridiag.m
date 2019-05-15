% Algorithm for solving tridiagonal linear systems of algebraic equations of 
% the form Ax=b via using Doolittle's LU factorization without selecting pivot.
% Here A is a tridiagonal matrix
% 10170437 Mark Taylor
function [L,U,x] = solveTridiag(A,b) 

% [L,U]=D_LU(A);    % Implementation of Doolittle's LU factorization can be
% replaced by following steps due to tridiagonal feature that A possesses.
 
% *****
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
elseif length(b)~=m
    error('The number of rows of A is not equal to length of column vector b!')
end

L=eye(n);
U=zeros(n);
U(1,1)=A(1,1);
for  i=2:n
   L(i,i-1)=A(i,i-1)/U(i-1,i-1);
   U(i,i)=A(i,i)-L(i,i-1)*A(i-1,i);
   U(i-1,i)=A(i-1,i);
end
% *****

% Solve Ly=b
y=zeros(n,1);
y(1)=b(1);
for i=2:n
    y(i)=b(i)-L(i,i-1)*y(i-1);
end

% Solve Ux=y
x=zeros(n,1);
x(n)=y(n)/U(n,n);
for i=n-1:-1:1
    x(i)=(y(i)-A(i,i+1)*x(i+1))/U(i,i);
    
end

% Only output x when there is one output argument
if nargout==1
    L=x;
end

end



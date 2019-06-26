% Cholesky Factorization, the decomposition of a symmetric positive definite  
% matrix into the product of a lower triangular matrix and its transpose.
% 10170437 Mark Taylor
function [x,L]= Cholesky(A,b)
    
[m,n]=size(A);
if m ~= n
    error('A must be symmetric positive definite!')
elseif isequal(A,A')==false || allPositive(eig(A))==false
    error('A must be symmetric positive definite!')
end

L=zeros(n);
% Perform Cholesky's method with no need to pivoting
% since A is symmetric positive definite.
for j=1:n
    L(j,j)=sqrt(A(j,j)-L(j,1:j-1)*L(j,1:j-1).');
    
    for i=j+1:n
        
        L(i,j)=(A(i,j)-L(i,1:j-1)*L(j,1:j-1).')/L(j,j);
        
    end   
    
end

if nargin==2 
    
% A=L*L.', L*y=b, y=L.'*x, => A*x=b
y=solveTril(L,b);
x=solveTriu(L.',y);

end

end
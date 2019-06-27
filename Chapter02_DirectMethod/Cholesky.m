% Cholesky Factorization, the decomposition of a symmetric positive definite  
% matrix into the product of a lower triangular matrix and its transpose.
% 10170437 Mark Taylor
function [L,x]= Cholesky(A,b)
    
[m,n]=size(A);
if m ~= n || isequal(A,A.')==false
    error('A must be symmetric positive definite!')
end

L=zeros(n);
% Perform Cholesky's method with no need to pivoting, since A is 
% symmetric positive definite, which determines its stability of
% eliminations in natural order.
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

if nargout<=1
    L=x;
end

end

end
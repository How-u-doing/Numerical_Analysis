% Cholesky factorization, the decomposition of a symmetric, positive definite  
% matrix into the product of a lower triangular matrix and its conjugate transpose.
% 10170437 Mark Taylor
function L= Cholesky(A)
    
[m,n]=size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
elseif isequal(A,A')==false || allPositive(eig(A))==false
    error('Input matrix must be symmetric & positive definite!')
end

L=zeros(n);
% Perform Cholesky's Method with no need to select
% pivot since A is symmetric, positive definite.
for j=1:n
    L(j,j)=sqrt(A(j,j)-L(j,1:j-1)*L(j,1:j-1).');
    
    for i=j+1:n
        
        L(i,j)=(A(i,j)-L(i,1:j-1)*L(j,1:j-1).')/L(j,j);
        
    end
    
    
end

end





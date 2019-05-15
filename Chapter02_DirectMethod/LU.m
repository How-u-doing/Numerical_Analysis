% LU matrix factorizaton, PA=LU, where L is a unit lower triangular matrix,
% U is an upper triangular matrix, and P is a permutation matrix.
% 10170437 Mark Taylor
function [L,U,P]=LU(A)% Require A is a nonsigular matrix
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input error! Input must be a square matrix!')
end

P=eye(n);
L=zeros(n);
for j=1:n-1 
    k=maxIndex(A(:,j),j,m);
    if A(k,j)~=0
        if k~=j
            % exchange j-th row and k-th row of identity matrix I
            I=eye(n);I(k,k)=0;I(k,j)=1;I(j,k)=1;I(j,j)=0;% construct P(j,k)
            P=I*P; %   P=P(n-1)...P(2)P(1)
            
            temp=A(j,:);     
            A(j,:)=A(k,:); % <-> A=P(k,j)*A
            A(k,:)=temp;
        end 
    else 
        error('Input is a singular matrix!')
    end    
    
% Perform Gaussian elimination and store the multiple coefficients l_ij
% in the lower triangular part of A so that the permutations of L(j+1:m,j)
% will be implemented as the pivot selecting process continues.
    L(j+1:m,j)=A(j+1:m,j)/A(j,j);
    A(j+1:m,j)=L(j+1:m,j);  % A(:,j) <= P(n-1)...P(j+1)L(:,j)
    A(j+1:m,j+1:n)=A(j+1:m,j+1:n)-L(j+1:m,j).*A(j,j+1:n);
    
end 

% Extract L and U from A
L = tril(A) - diag(diag(A)) + eye(n);  % L=tril(A,-1)+eye(n)
U = triu(A);
        
end



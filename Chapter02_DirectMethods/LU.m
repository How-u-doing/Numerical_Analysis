% LU Factorization, PA=LU, where L is a unit lower triangular matrix,
% U is an upper triangular matrix, and P is a permutation matrix.
% 10170437 Mark Taylor

function [L,U,P]=LU(A)                                   % Require A is invertible

[m,n] = size(A);
if m ~= n
    error('Input must be square!')
end

P=eye(n);
L=zeros(n);
for j=1:n-1 
    k=maxIndex(A(:,j),j,n);
    if abs(A(k,j))>eps
        if k~=j
            % Exchange j-th row with k-th row of identity matrix I
            I=eye(n);I(k,k)=0;I(k,j)=1;I(j,k)=1;I(j,j)=0;% Construct P(j,k)
            P=I*P;                                       % P=P(n-1)...P(2)P(1)
            
            % Note the permutations of l_ijs.
            temp=A(j,:);
            A(j,:)=A(k,:);
            A(k,:)=temp;
        end 
    else 
        error('A is singular!')
    end    
    
    % Perform Gaussian elimination and store the multipliers i.e. l_ijs in
    % the lower triangular part of A so that the permutations of L(j+1:n,j)
    % will be implemented as the partial pivoting process continues.
    L(j+1:n,j)=A(j+1:n,j)/A(j,j);
    A(j+1:n,j)=L(j+1:n,j);                               % A(:,j) <= P(n-1)...P(j+1)L(:,j)
    A(j+1:n,j+1:n)=A(j+1:n,j+1:n)-L(j+1:n,j).*A(j,j+1:n);
    
end 

% Extract L and U from A
L = tril(A) - diag(diag(A)) + eye(n);                    % L=tril(A,-1)+eye(n)
U = triu(A);
        
end



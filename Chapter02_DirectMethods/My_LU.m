% LU Factorization, PA=LU, where L is a unit lower triangular matrix,
% U is an upper triangular matrix, and P is a permutation matrix.
% 10170437 Mark Taylor

function [L,U,P]=My_LU(A)                                % Require A is invertible

[m,n] = size(A);
if m ~= n
    error('Input must be square!')
end

P=eye(n);
L=eye(n);
U=A;
K=zeros(1,n-1);
for j=1:n-1
    k=maxIndex(U(:,j),j,n);
    K(j)=k;
    if abs(U(k,j))>eps 
        if k~=j
            % Exchange j-th row with k-th row of identity matrix I
            I=eye(n);I(k,k)=0;I(k,j)=1;I(j,k)=1;I(j,j)=0;% Construct P(j,k)
            P=I*P;                                       % P=P(n-1)...P(2)P(1)
            
            temp=U(j,j:n);
            U(j,j:n)=U(k,j:n);                           % <-> U=P(k,j)*U or U([j,k],j:n)=U([k,j],j:n)
            U(k,j:n)=temp;
            
% Now let's analyze the efficiency of U([j,k],j:n)=U([k,j],j:n).
% In my opinion,calling or using U([k,j],j:n) to perform the assignment may
% construct two copies of row vector (or a copy of 2-by-(n-j+1) matrix) to 
% avoid chaos, which might cause more consumption of memory especially when 
% n is extremely large.
% Of course, assignment should be optimized by the compiler.


        end 
    else 
        error('A is singular!')
    end    
    
    % Perform Gaussian elimination
    L(j+1:n,j)=U(j+1:n,j)/U(j,j);
    U(j+1:n,j:n)=U(j+1:n,j:n)-L(j+1:n,j).*U(j,j:n);    
    % /******  <-> following part
    % 
    % ej=zeros(1,n); ej(j)=1;                            % unit row vector e(j)
    % L(j+1:n,j)=U(j+1:n,j)/U(j,j);
    % G=eye(n)-L(:,j)*ej;    
    % U=G*U;
    %
    % ******/ 
    
end

%                                                                                  _
% Until here, the construction of L is not finished yet owing to the fact that inv(Gk)=
% P(n-1)...P(k+1)inv(Gk)P(k+1)...P(n-1).Therefore,L(:,k) <= P(n-1)...P(k+1)L(:,k) 
S=P;
for j=1:n-2
    t=K(j);
    if t~=j
        I=eye(n);I(t,t)=0;I(t,j)=1;I(j,t)=1;I(j,j)=0;    % construct P(j,t)
        S=S*I;                                           % S=P(n-1)...P(2)P(1) * P(1)P(2)...P(n-2) 
    end
    
    L(:,j)=S*L(:,j);                                     % permutate L(:,j) or rather L(j+1:m,j)
end

% Indeed, another alternative(should be also a better & practical solution)
% is to store L(j+1:n,j) in A for each j from 1 to n-1. Then permutations 
% of L(j+1:m,j) will be implemented as the pivot selecting process continues.
% Ultimately, what we should do is to extract L & U respectively from the
% lower and upper triangular part of A.  
% For example, L(:,1)=P(n-1)...P(2)L(:,1). 
% See "LU.m" for further details.


end


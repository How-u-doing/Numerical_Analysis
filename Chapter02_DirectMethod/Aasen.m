% Aasen's algorithm :  PAP'=LTL',  where A is a n-by-n (dense) symmetric
% matrix, L is a unit lower triangular matrix,P is a permutation matrix, 
% and T is a symmetric tridiagonal matrix.
% 10170437 Mark Taylor
function [P,L,T]=Aasen(A)

[m,n]=size(A);
if m ~= n || isequal(A,A')==false 
    error('Input error! Input must be a symmetric matrix!')
end

P=eye(n);
L=zeros(n);
for i=1:n-2
    r=maxIndex(A(:,i),i+1,n); 
    if abs(A(r,i))>eps
        if r~=i+1
            A([i+1,r],:)=A([r,i+1],:); % row (i+1) <-> row r
            A(:,[i+1,r])=A(:,[r,i+1]); % column (i+1) <-> column r
            I=eye(n);I(r,r)=0;I(i+1,r)=1;I(r,i+1)=1;I(i+1,i+1)=0;% construct P(i+1,r)
            P=I*P; %   P=P(n-2)...P(2)P(1)
        end
    
    end
    
% Perform Gaussian elimination and store the multiple coefficients l_ij
% in the lower & upper triangular part of A so that the permutations of 
% L(r+2:n,i) & L(r+2:n,i).'  will be implemented as the pivot selecting 
% process continues.

%      _ _
% A <= GAG' ... until A is tridiagonal 
    L(i+2:n,i)=A(i+2:n,i)/A(i+1,i);
    A(i+2:n,i)=L(i+2:n,i);  % {A(:,i)}(i+2:n) <= {P(n-2)...P(i+1)L(:,i)}(i+2:n)
    A(i+2:n,i+1:n)=A(i+2:n,i+1:n)-L(i+2:n,i).*A(i+1,i+1:n);    
    
    % Transpose above steps
    L(i,i+2:n)=A(i,i+2:n)/A(i,i+1);
    A(i,i+2:n)=L(i,i+2:n);  % {A(i:,)}(i+2:n) <= {L(:,i)P(i+1)...P(n-2)}(i+2:n)
    A(i+1:n,i+2:n)=A(i+1:n,i+2:n)-A(i+1:n,i+1).*L(i,i+2:n);    
    
end

% Extract L and T from A

L=zeros(n);
% Shift this sub-lower-triangular part(i>=j+2) right a single unit
L(3:n,2:n-1)=tril(A(3:n,1:n-2));
L=L+eye(n);

T=A-tril(A,-2)-triu(A,2);

end



% Givens Rotation &  QR factorization
% 10170437 Mark Taylor
function [Q,R]=G_QR(A)
[m,n] = size(A);
if m ~= n % make sure A is a square matrix to continue following steps
    error('Input must be a square matrix!')
end
copy_A=A;
Q=eye(n);
for j=1:n-1
    for i=j+1:n
        if A(i,j)~=0
            T=Givens(A(j,j),A(i,j),j,i,n);
            A=T*A;
            Q=T*Q;% Q=T(n,n-1)T(n-2,n)...T(1,n)...T(1,2)
        else 
            continue;
        end
    end
end
R=Q*copy_A;
Q=Q.';
end

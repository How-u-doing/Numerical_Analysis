% Givens Rotation &  QR factorization
% 10170437 Mark Taylor
function [Q,R]=G_QR(A)
[m,n] = size(A);
if m ~= n 
    error('Input must be a square matrix!')
end
R=A;
Q=eye(n);
for j=1:n-1
    for i=j+1:n
        if R(i,j)~=0
            T=Givens(R(j,j),R(i,j),j,i,n);
            R=T*R;     % R=T(n,n-1)T(n-2,n)...T(1,n)...T(1,2)A
            Q=T*Q;     % Q=T(n,n-1)T(n-2,n)...T(1,n)...T(1,2)
        else 
            continue;
        end
    end
end

Q=Q.';
end

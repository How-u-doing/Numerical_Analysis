% Givens QR Factorization
% 10170437 Mark Taylor
function [Q,R]=G_QR(A)
[m,n] = size(A);
if m < n 
    error('The number of rows must be no less than the number of columns!')
elseif m==n
    n=n-1;
end
R=A;
Q=eye(m);

for j=1:n
    for i=j+1:m
        if R(i,j)~=0
            T=GV(R(j,j),R(i,j),j,i,m);
            R=T*R;     % m=n: R=T(n-1,m)T(n-2,m)...T(1,m)...T(1,2)A
            Q=T*Q;     % m=n: Q=T(1,2)...T(1,m)...T(n-2,m)T(n-1,m)
        else 
            continue;
        end
    end
end

Q=Q.';
end

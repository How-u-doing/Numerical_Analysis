% QR decomposition via Givens Rotation 
% 10170437 Mark Taylor
function [Q,R]=givens_QR(A)
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
            G=givens_cs(R(j,j),R(i,j));
            R([j,i],:)=G*R([j,i],:);
            Q([j,i],:)=G*Q([j,i],:);
            
        else 
            continue;
        end
    end
end

Q=Q.';
end
